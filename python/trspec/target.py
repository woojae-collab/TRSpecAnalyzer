"""Target analysis: fit a kinetic scheme (compartment model) of your own.

Port of parseTargetScheme, targetSchemeText, targetRateMatrix, targetConc
and fitTargetAnalysis of TRSpecAnalyzer.m.

Indices are 0-based here: species 0 .. n_s-1, and a step that ends in the
ground state has  to == -1.
"""
import re
import numpy as np
from scipy.optimize import nnls
from ._optim import fminsearch
from ._util import Result, lstsq_min_norm
from .kinetics import (ca_basis, ca_options, exp_irf_conv, irf_delay_step,
                       irf_unresolved, irf_width_lost)
from .uncertainty import fit_param_uncertainty

_SINK_NAMES = ('gs', 's0', '0', 'ground')
_NUM = r'([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)'
_NAM = r'([A-Za-z][A-Za-z0-9_]*)'
_REALMIN = np.finfo(float).tiny


def _str2double(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float('nan')


def _is_sink(name):
    return name.lower() in _SINK_NAMES


def _read_range(text):
    r2 = text.replace('..', ' ')
    r2 = re.sub(r'(?<=[0-9.])\s*-\s*(?=[0-9.])', ' ', r2)
    r2 = r2.replace(',', ' ')
    vals = []
    for tok in r2.split():
        try:
            vals.append(float(tok))
        except ValueError:
            break
    if any(not np.isfinite(v) for v in vals):
        return []
    return vals


def parse_target_scheme(text):
    """Read a kinetic scheme written as text.

    One statement per line (or separated by ';'); '%' and '#' start a
    comment. Examples:
        A -> B                  step with its own rate (k1, k2, ... by order)
        A -> B : k1             step with a named rate
        B -> C : 2*k1           rate tied to another one by a factor
        A -> B -> C             chain
        A <-> B : kf, kb        reversible step
        B -> 2 C                two products per event
        A -> GS                 loss to the ground state (GS, S0, 0, ground)
        species: A, B, C        order of the compartments
        init: A = 0.7, B = 0.3  initial populations (default: first one = 1)
        zero: B 700 900         spectrum of B is zero from 700 to 900
        equal: C = B            spectra equal (optionally a factor and a range:
        equal: C = 0.5*B 600 700)
        nonneg                  all spectra >= 0   (nonneg: A, B  for some)

    Returns a Result with species, n_s, params, n_p, steps (list of dicts:
    from, to, coef, par, mult, label, expr), c0, zero, rel, nonneg, sink_name.
    """
    if not isinstance(text, str):
        text = '\n'.join(str(c) for c in text)
    stm = re.split(r'[\n\r;]+', text)
    state = {'sink_name': '', 'n_auto': 0}
    species, params, steps = [], [], []
    init_txt, zero_txt, equal_txt, nonneg_txt = [], [], [], []
    nonneg_all = False

    def check_name(name, line):
        if not name or not re.match(r"^[A-Za-z0-9_()'^{}]+$", name):
            raise ValueError('Scheme line "%s": "%s" is not a valid compartment name '
                             "(letters, digits, _ ( ) ' only)." % (line, name))
        if not _is_sink(name) and re.match(r'^[0-9.]+$', name):
            raise ValueError('Scheme line "%s": a compartment name cannot be a number ("%s").'
                             % (line, name))

    def species_index(name):
        if name not in species:
            species.append(name)
        return species.index(name)

    def parse_term(term, line):
        term = term.strip()
        coef = 1.0
        tk = re.match(r'^([0-9]*\.?[0-9]+)\s*\*\s*(\S+)$', term)
        if not tk:
            tk = re.match(r'^([0-9]*\.?[0-9]+)\s+(\S+)$', term)
        if tk:
            coef = _str2double(tk.group(1))
            name = tk.group(2)
        else:
            name = term
        check_name(name, line)
        if not np.isfinite(coef) or coef <= 0:
            raise ValueError('Scheme line "%s": bad factor in "%s".' % (line, term))
        return name, coef

    def parse_rate(s, line):
        s = s.strip()
        mult = 1.0
        t1 = re.match('^' + _NAM + '$', s)
        t2 = re.match('^' + _NUM + r'\s*\*\s*' + _NAM + '$', s)
        t3 = re.match('^' + _NAM + r'\s*\*\s*' + _NUM + '$', s)
        t4 = re.match('^' + _NAM + r'\s*/\s*' + _NUM + '$', s)
        if t1:
            pname = t1.group(1)
        elif t2:
            pname, mult = t2.group(2), _str2double(t2.group(1))
        elif t3:
            pname, mult = t3.group(1), _str2double(t3.group(2))
        elif t4:
            den = _str2double(t4.group(2))
            pname, mult = t4.group(1), (1.0 / den if den != 0 else float('inf'))
        else:
            raise ValueError('Scheme line "%s": cannot read the rate "%s" (use  k1,  2*k1  or  k1/2).'
                             % (line, s))
        if not np.isfinite(mult) or mult <= 0:
            raise ValueError('Scheme line "%s": the factor in "%s" must be positive.' % (line, s))
        return pname, mult

    def add_step(src, dst, coef, rate_txt, line):
        i_from = species_index(src)
        if _is_sink(dst):
            i_to = -1
            if not state['sink_name']:
                state['sink_name'] = dst
        else:
            i_to = species_index(dst)
        for st in steps:
            if st['from'] == i_from and st['to'] == i_to:
                raise ValueError('The step "%s -> %s" appears twice.' % (src, dst))
        if not rate_txt:
            state['n_auto'] += 1
            pname, mlt = '#%d' % state['n_auto'], 1.0
        else:
            pname, mlt = parse_rate(rate_txt, line)
        if pname not in params:
            params.append(pname)
        ip = params.index(pname)
        lab = '%s -> %s' % (src, dst) if coef == 1 else '%s -> %g %s' % (src, coef, dst)
        steps.append({'from': i_from, 'to': i_to, 'coef': coef, 'par': ip,
                      'mult': mlt, 'label': lab, 'expr': ''})

    for raw in stm:
        s = re.sub(r'[%#].*$', '', raw).strip()
        if not s:
            continue
        tk_key = re.match(r'(?i)^(init|zero|equal|nonneg|species)\s*:\s*(.*)$', s)
        if tk_key:
            arg = tk_key.group(2).strip()
            key = tk_key.group(1).lower()
            if key == 'init':
                init_txt.append(arg)
            elif key == 'zero':
                zero_txt.append(arg)
            elif key == 'equal':
                equal_txt.append(arg)
            elif key == 'nonneg':
                if not arg:
                    nonneg_all = True
                else:
                    nonneg_txt.append(arg)
            else:
                for nm in re.split(r'\s*,\s*', arg):
                    if not nm:
                        continue
                    check_name(nm, s)
                    if _is_sink(nm):
                        raise ValueError('Scheme line "%s": the ground state (%s) is not a compartment.'
                                         % (s, nm))
                    species_index(nm)
            continue
        if s.lower() == 'nonneg':
            nonneg_all = True
            continue
        rate_lab = []
        rx = s
        ic = s.find(':')
        if ic >= 0:
            rx = s[:ic].strip()
            rate_lab = [x for x in re.split(r'\s*,\s*', s[ic + 1:].strip()) if x]
        parts = re.split(r'<->|->', rx)
        arrows = re.findall(r'<->|->', rx)
        if not arrows:
            raise ValueError('Scheme line "%s": expected "A -> B".' % s)
        n_here = len(arrows) + sum(1 for a in arrows if a == '<->')
        if rate_lab and len(rate_lab) != n_here:
            raise ValueError('Scheme line "%s": %d step(s) but %d rate label(s) after ":".'
                             % (s, n_here, len(rate_lab)))
        names, coefs = [], []
        for part in parts:
            nm, cf = parse_term(part, s)
            names.append(nm)
            coefs.append(cf)
        i_lab = 0
        for aa, arrow in enumerate(arrows):
            src, dst = names[aa], names[aa + 1]
            if _is_sink(src):
                raise ValueError('Scheme line "%s": the ground state (%s) cannot be the source of a step.'
                                 % (s, src))
            if src == dst:
                raise ValueError('Scheme line "%s": a step needs two different compartments.' % s)
            i_lab += 1
            add_step(src, dst, coefs[aa + 1], rate_lab[i_lab - 1] if rate_lab else '', s)
            if arrow == '<->':
                if _is_sink(dst):
                    raise ValueError('Scheme line "%s": a reversible step cannot involve the ground state.' % s)
                if coefs[aa + 1] != 1:
                    raise ValueError('Scheme line "%s": no stoichiometric factor on a reversible step.' % s)
                i_lab += 1
                add_step(dst, src, 1.0, rate_lab[i_lab - 1] if rate_lab else '', s)

    n_s = len(species)
    if n_s < 1:
        raise ValueError('The scheme has no compartments. Example:  A -> B')
    nxt = 1
    for pp in range(len(params)):
        if params[pp][0] != '#':
            continue
        while ('k%d' % nxt) in params:
            nxt += 1
        params[pp] = 'k%d' % nxt
    for st in steps:
        st['expr'] = params[st['par']] if st['mult'] == 1 else '%.10g*%s' % (st['mult'], params[st['par']])

    c0 = np.zeros(n_s)
    if not init_txt:
        c0[0] = 1.0
    else:
        for it in re.split(r'\s*,\s*', ','.join(init_txt)):
            it = it.strip()
            if not it:
                continue
            tk = re.match(r'^(\S+?)\s*=\s*(\S+)$', it)
            if not tk:
                raise ValueError('init: cannot read "%s" (use  init: A = 1).' % it)
            if tk.group(1) not in species:
                raise ValueError('init: "%s" is not a compartment of the scheme.' % tk.group(1))
            val = _str2double(tk.group(2))
            if not np.isfinite(val) or val < 0:
                raise ValueError('init: the population of %s must be >= 0.' % tk.group(1))
            c0[species.index(tk.group(1))] = val
        if not np.any(c0 > 0):
            raise ValueError('init: at least one compartment needs a population > 0.')

    zero = []
    for z in zero_txt:
        tk = re.match(r'^(\S+)\s+(.*)$', z)
        if not tk:
            raise ValueError('zero: cannot read "%s" (use  zero: B 700 900).' % z)
        if tk.group(1) not in species:
            raise ValueError('zero: "%s" is not a compartment of the scheme.' % tk.group(1))
        v = _read_range(tk.group(2))
        if len(v) != 2:
            raise ValueError('zero: "%s" needs two numbers (from, to).' % z)
        zero.append({'sp': species.index(tk.group(1)), 'lo': min(v), 'hi': max(v)})

    rel = []
    for q in equal_txt:
        ie = q.find('=')
        tk = None
        if ie >= 0:
            rest = re.sub(r'\s*\*\s*', '*', q[ie + 1:].strip())
            tk = re.match(r'^(\S+)\s*(.*)$', rest)
        if not tk:
            raise ValueError('equal: cannot read "%s" (use  equal: C = B  or  equal: C = 0.5*B 600 700).' % q)
        dep_n = q[:ie].strip()
        mas_n = tk.group(1)
        rt = 1.0
        tm = re.match('^' + _NUM + r'\*(.+)$', mas_n)
        if tm:
            rt = _str2double(tm.group(1))
            mas_n = tm.group(2)
        if dep_n not in species or mas_n not in species:
            raise ValueError('equal: "%s" names a compartment that is not in the scheme.' % q)
        i_d, i_m = species.index(dep_n), species.index(mas_n)
        if i_d == i_m:
            raise ValueError('equal: "%s" relates a compartment to itself.' % q)
        if not np.isfinite(rt) or rt == 0:
            raise ValueError('equal: the factor in "%s" must be a non-zero number.' % q)
        lo, hi = -np.inf, np.inf
        if tk.group(2).strip():
            v = _read_range(tk.group(2))
            if len(v) != 2:
                raise ValueError('equal: the range in "%s" needs two numbers.' % q)
            lo, hi = min(v), max(v)
        for r in rel:
            if r['dep'] == i_d and max(r['lo'], lo) <= min(r['hi'], hi):
                raise ValueError('equal: "%s" is tied twice in overlapping ranges.' % species[i_d])
        rel.append({'dep': i_d, 'master': i_m, 'ratio': rt, 'lo': lo, 'hi': hi})
    for r in rel:
        cur = r['master']
        for _ in range(n_s):
            nx = next((x for x in rel if x['dep'] == cur), None)
            if nx is None:
                break
            cur = nx['master']
            if cur == r['dep']:
                raise ValueError('equal: the constraints on "%s" are circular.' % species[r['dep']])

    nonneg = np.zeros(n_s, dtype=bool)
    if nonneg_all:
        nonneg[:] = True
    for item in nonneg_txt:
        for nm in re.split(r'\s*,\s*', item):
            if not nm:
                continue
            if nm not in species:
                raise ValueError('nonneg: "%s" is not a compartment of the scheme.' % nm)
            nonneg[species.index(nm)] = True

    return Result(species=species, n_s=n_s, params=params, n_p=len(params), steps=steps,
                  c0=c0, zero=zero, rel=rel, nonneg=nonneg,
                  sink_name=state['sink_name'] or 'GS')


def target_scheme_text(sch):
    """The scheme written out as text lines (canonical form)."""
    sp = sch.species
    lines = ['species: ' + ', '.join(sp)]
    for st in sch.steps:
        dst = sch.sink_name if st['to'] < 0 else sp[st['to']]
        if st['coef'] != 1:
            dst = '%g %s' % (st['coef'], dst)
        lines.append('%s -> %s : %s' % (sp[st['from']], dst, st['expr']))
    c0 = np.asarray(sch.c0, dtype=float).ravel()
    default = np.zeros(len(sp))
    default[0] = 1
    if not np.array_equal(c0, default):
        items = ['%s = %.10g' % (sp[j], c0[j]) for j in range(len(sp)) if c0[j] != 0]
        lines.append('init: ' + ', '.join(items))
    for z in sch.zero:
        lines.append('zero: %s %.10g %.10g' % (sp[z['sp']], z['lo'], z['hi']))
    for r in sch.rel:
        rhs = sp[r['master']] if r['ratio'] == 1 else '%.10g*%s' % (r['ratio'], sp[r['master']])
        ln = 'equal: %s = %s' % (sp[r['dep']], rhs)
        if np.isfinite(r['lo']) and np.isfinite(r['hi']):
            ln = '%s %.10g %.10g' % (ln, r['lo'], r['hi'])
        lines.append(ln)
    nonneg = np.asarray(sch.nonneg, dtype=bool)
    if nonneg.any():
        if nonneg.all():
            lines.append('nonneg')
        else:
            lines.append('nonneg: ' + ', '.join(s for s, f in zip(sp, nonneg) if f))
    return lines


def target_rate_matrix(sch, tau_par):
    """Rate matrix K (dc/dt = K c) for the time constant of every rate
    parameter; a step runs with  mult / tau."""
    tau_par = np.asarray(tau_par, dtype=float).ravel()
    K = np.zeros((sch.n_s, sch.n_s))
    for st in sch.steps:
        kk = st['mult'] / tau_par[st['par']]
        a, b = st['from'], st['to']
        K[a, a] -= kk
        if b >= 0:
            K[b, a] += st['coef'] * kk
    return K


def _rcond(v):
    try:
        c = np.linalg.cond(v, 1)
    except Exception:
        return 0.0
    return 0.0 if not np.isfinite(c) or c == 0 else 1.0 / c


def target_conc(t, K, c0, t0, fwhm, irf_mode='irf'):
    """Populations c(t) of the scheme (eigen-decomposition of K, every mode
    convolved with the Gaussian IRF).

    Returns (C, obs_tau, mode_w): C (n_t x n_s); the observable lifetimes
    (ascending, inf for a non-decaying mode) and mode_w with
    C = E * mode_w.T, E being the mode time profiles.
    """
    t = np.asarray(t, dtype=float).ravel()
    K = np.asarray(K, dtype=float)
    n = K.shape[0]
    c0 = np.asarray(c0, dtype=float).ravel()
    sc = np.max(np.abs(np.diag(K)))
    if not (sc > 0):
        sc = 1.0
    lam = V = None
    for attempt in range(5):
        Kp = K if attempt == 0 else K * (1 + 10.0 ** (attempt - 7) * np.arange(1, n + 1))[None, :]
        lam, V = np.linalg.eig(Kp)
        if _rcond(V) > 1e-9:
            break
    if np.any(np.abs(np.imag(lam)) > 1e-9 * sc):
        raise ValueError('The scheme gives oscillating populations (complex decay rates, e.g. a '
                         'closed loop of three or more compartments). This is not supported.')
    lam = np.real(lam)
    V = np.real(V)
    if np.any(lam > 1e-9 * sc):
        raise ValueError('The scheme makes the populations grow without bound '
                         '(a loop that contains a multiplying step?).')
    if _rcond(V) < 1e-12:
        raise ValueError('The rate matrix is degenerate; change the time constants slightly.')
    a = np.linalg.solve(V, c0)
    no_irf = isinstance(irf_mode, str) and irf_mode.lower() == 'none'
    B = np.zeros((t.size, n))
    tau_m = np.full(n, np.inf)
    for m in range(n):
        if lam[m] < -1e-12 * sc:
            tau_m[m] = -1.0 / lam[m]
        if no_irf:
            x = t - t0
            b = np.zeros(x.shape)
            on = x >= 0
            b[on] = np.exp(-x[on] / tau_m[m]) if np.isfinite(tau_m[m]) else 1.0
        else:
            b = exp_irf_conv(t, tau_m[m], t0, fwhm)
        B[:, m] = b
    order = np.argsort(tau_m, kind='stable')
    mode_w = (V * a[None, :])[:, order]
    C = B[:, order] @ mode_w.T
    C[~np.isfinite(C)] = 0
    return C, tau_m[order], mode_w


def _nn_solve(Bm, y, nnc):
    if np.all(nnc):
        return nnls(Bm, y)[0]
    n1 = int(np.count_nonzero(nnc))
    n2 = int(np.count_nonzero(~nnc))
    zz = nnls(np.hstack([Bm[:, nnc], Bm[:, ~nnc], -Bm[:, ~nnc]]), y)[0]
    xa = np.zeros(nnc.size)
    xa[nnc] = zz[:n1]
    xa[~nnc] = zz[n1:n1 + n2] - zz[n1 + n2:]
    return xa


def fit_target_analysis(D, wl, t, sch, tau_init, tau_fixed=False, t0_init=0.0,
                        fwhm_init=0.15, t0_fixed=True, fwhm_fixed=True,
                        irf_mode='irf', t_range=None, tau_lo=None, tau_hi=None,
                        ca=None, fw_mode='auto'):
    """Fit D(lambda, t) = sum_s SAS_s(lambda) * c_s(t) for the scheme sch.

    sch       : scheme from parse_target_scheme (or its text)
    tau_init  : time constant of every rate parameter (order of sch.params)
    tau_fixed : bool per parameter (or one for all)
    irf_mode  : 'irf' (Gaussian IRF) or 'none'
    t_range   : optional (tmin, tmax) fit window
    tau_lo/hi : optional limits per parameter (NaN = none)
    ca        : coherent-artifact terms at time zero (see ca_options)
    fw_mode   : 'auto' | 'plain' | 'limit', handling of a fitted IRF FWHM
                (see fit_global_analysis)

    Returns a Result with tau, t0, fwhm, sas (NaN where not determinable),
    fit, rss, rms, iters, converged, initial_rms, n_used, tau_std, t0_std,
    fwhm_std (inf = not determined by the data), at_limit, conc, K, obs_tau,
    n_obs, das, n_free, t_sel, reordered, scheme, ca_order, ca_spec,
    ca_fwhm, ca_fwhm_std, ca_mode, fwhm_lo, fwhm_at_limit, irf_retry.
    """
    fw_mode = fw_mode or 'auto'
    if isinstance(sch, (str, list, tuple)):
        sch = parse_target_scheme(sch)
    D = np.asarray(D, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    wl = np.asarray(wl, dtype=float).ravel()
    M, N = D.shape
    n_s, n_p = sch.n_s, sch.n_p
    tau_init = np.atleast_1d(np.asarray(tau_init, dtype=float)).ravel()
    tau_fixed = np.atleast_1d(np.asarray(tau_fixed, dtype=bool)).ravel()
    if tau_fixed.size == 1:
        tau_fixed = np.repeat(tau_fixed, n_p)
    if tau_init.size != n_p or tau_fixed.size != n_p:
        raise ValueError('fit_target_analysis: one time constant per rate parameter (%d) is needed.' % n_p)
    if np.any(~np.isfinite(tau_init)) or np.any(tau_init <= 0):
        raise ValueError('All time constants must be positive and finite.')
    if t.size != N or wl.size != M:
        raise ValueError('fit_target_analysis: axes do not match the data size.')
    tau_lo = np.full(n_p, np.nan) if tau_lo is None or np.size(tau_lo) == 0 else np.asarray(tau_lo, dtype=float).ravel()
    tau_hi = np.full(n_p, np.nan) if tau_hi is None or np.size(tau_hi) == 0 else np.asarray(tau_hi, dtype=float).ravel()
    if tau_lo.size != n_p or tau_hi.size != n_p:
        raise ValueError('fit_target_analysis: limits need one entry per rate parameter.')
    no_irf = isinstance(irf_mode, str) and irf_mode.lower() == 'none'
    if no_irf:
        t0_fixed = True
        fwhm_fixed = True
    elif not (fwhm_init > 0):
        raise ValueError('The IRF FWHM must be positive.')
    t0_fixed_arg, fwhm_fixed_arg = bool(t0_fixed), bool(fwhm_fixed)
    ca = ca_options(ca)
    if no_irf:
        ca['order'] = 0
    ca_k = ca['order']
    ca_own = ca_k > 0 and ca['mode'] != 'irf'
    ca_free = ca_k > 0 and ca['mode'] == 'free'
    ca_w0 = ca['fwhm']
    if ca_own and not (ca_w0 > 0):
        ca_w0 = float(fwhm_init)
    t_sel = np.ones(N, dtype=bool)
    if t_range is not None and np.size(t_range) == 2 and np.all(np.isfinite(t_range)):
        t_sel = (t >= np.min(t_range)) & (t <= np.max(t_range))
    if not t_sel.any():
        raise ValueError('No delay points inside the fit t range.')
    t_span = t[t_sel].max() - t[t_sel].min()
    if not (t_span > 0):
        t_span = np.max(np.abs(t)) + np.finfo(float).eps
    fw_res, fw_lo, fw_wall0 = 0.0, 0.0, False
    if not fwhm_fixed:
        dz_irf = irf_delay_step(t[t_sel], t0_init)
        if np.isfinite(dz_irf) and dz_irf > 0:
            fw_res = 0.5 * dz_irf
            fw_lo = min(fw_res, fwhm_init) * (1 - 1e-9)
            fw_wall0 = fw_mode == 'limit'
    if ca_k > 0:
        w0 = ca_w0 if ca_own else fwhm_init
        n_near = int(np.count_nonzero(t_sel & (np.abs(t - t0_init) <= 1.5 * w0)))
        if n_near < ca_k + 2:
            raise ValueError('Coherent artifact: %d fitted delay(s) lie within 1.5 FWHM of '
                             'time zero, %d are needed for %d term(s). Use fewer terms or '
                             'switch them off.' % (n_near, ca_k + 2, ca_k))

    # spectral constraints -> groups of wavelengths with the same pattern
    code = np.zeros((M, n_s), dtype=int)
    for zz, r in enumerate(sch.rel):
        rows = (wl >= r['lo']) & (wl <= r['hi'])
        code[rows, r['dep']] = zz + 1
    for z in sch.zero:
        rows = (wl >= z['lo']) & (wl <= z['hi'])
        code[rows, z['sp']] = -1
    pat, grp = np.unique(code, axis=0, return_inverse=True)
    grp = np.asarray(grp).ravel()
    n_g = pat.shape[0]
    Tg, nnG = [], []
    nonneg = np.asarray(sch.nonneg, dtype=bool)
    for g in range(n_g):
        Tm = np.eye(n_s)
        Tm[pat[g, :] == -1, :] = 0
        for _ in range(n_s):
            for s in np.flatnonzero(pat[g, :] > 0):
                rr = sch.rel[pat[g, s] - 1]
                Tm[s, :] = rr['ratio'] * Tm[rr['master'], :]
        keep_c = np.any(Tm != 0, axis=0)
        Tg.append(Tm[:, keep_c])
        nnG.append(nonneg[keep_c])
    row_finite = np.all(np.isfinite(D[:, t_sel]), axis=1)

    idx_free = np.flatnonzero(~tau_fixed)
    n_free = idx_free.size
    g_lo, g_hi = 1e-7 * t_span, 1e3 * t_span
    has_lim = np.zeros(n_free, dtype=bool)
    l_lo = np.zeros(n_free)
    l_hi = np.zeros(n_free)
    for q in range(n_free):
        lo, hi = tau_lo[idx_free[q]], tau_hi[idx_free[q]]
        has_lim[q] = (np.isfinite(lo) and lo > 0) or (np.isfinite(hi) and hi > 0)
        if not (np.isfinite(lo) and lo > 0):
            lo = g_lo
        if not (np.isfinite(hi) and hi > 0):
            hi = g_hi
        if not (hi > lo):
            raise ValueError('Parameter %s: the upper limit must exceed the lower limit.'
                             % sch.params[idx_free[q]])
        l_lo[q], l_hi[q] = np.log(lo), np.log(hi)
        if has_lim[q] and (tau_init[idx_free[q]] < lo or tau_init[idx_free[q]] > hi):
            raise ValueError('Parameter %s: the starting value %g is outside its limits [%g, %g].'
                             % (sch.params[idx_free[q]], tau_init[idx_free[q]], lo, hi))

    def tau2x(tv):
        xv = np.zeros(n_free)
        for q1 in range(n_free):
            lt = np.log(max(tv[q1], _REALMIN))
            if has_lim[q1]:
                fr = (lt - l_lo[q1]) / (l_hi[q1] - l_lo[q1])
                fr = min(max(fr, 0.0), 1.0)
                xv[q1] = np.arcsin(2 * fr - 1)
            else:
                xv[q1] = lt
        return xv

    x0 = list(tau2x(tau_init[idx_free]))
    if not t0_fixed:
        x0.append(float(t0_init))
    if not fwhm_fixed:
        x0.append(np.log(fwhm_init))
    if ca_free:
        x0.append(np.log(ca_w0))
    x0 = np.asarray(x0, dtype=float)

    st = {'tau': tau_init.copy(), 't0': float(t0_init), 'fwhm': float(fwhm_init),
          'C': np.zeros((N, n_s)), 'n_loss': 0, 'n_lin': 0,
          'row_done': np.zeros(M, dtype=bool), 'ca_w': ca_w0,
          'Cca': np.zeros((N, ca_k)), 'Aca': np.zeros((M, ca_k)),
          'fw_wall': fw_wall0}
    t_min, t_max = t.min(), t.max()

    def unpack(x):
        tau_cur = tau_init.copy()
        for q2 in range(n_free):
            if has_lim[q2]:
                lt = l_lo[q2] + (l_hi[q2] - l_lo[q2]) * (np.sin(x[q2]) + 1) / 2
            else:
                lt = x[q2]
            with np.errstate(over='ignore'):
                tau_cur[idx_free[q2]] = np.exp(lt)
        st['tau'] = tau_cur
        p = n_free
        if not t0_fixed:
            st['t0'] = float(x[p])
            p += 1
        if not fwhm_fixed:
            with np.errstate(over='ignore'):
                st['fwhm'] = float(np.exp(x[p]))
            p += 1
        if ca_free:
            with np.errstate(over='ignore'):
                st['ca_w'] = float(np.exp(x[p]))

    def objective(x):
        x = np.asarray(x, dtype=float).ravel()
        unpack(x)
        A = np.zeros((M, n_s))
        fit_m = np.zeros((M, N))
        loss = 1e30
        st['n_lin'] = 0
        row_done = np.zeros(M, dtype=bool)
        st['row_done'] = row_done
        Aca = np.zeros((M, ca_k))
        st['Aca'] = Aca
        tau_cur, t0_cur, fwhm_cur = st['tau'], st['t0'], st['fwhm']
        ca_w_cur = st['ca_w']
        if n_free > 0:
            tf = tau_cur[idx_free]
            if np.any(~np.isfinite(tf)) or np.any(tf > g_hi * (1 + 1e-9)) or np.any(tf < g_lo * (1 - 1e-9)):
                return loss, A, fit_m
        if not no_irf:
            if (not np.isfinite(fwhm_cur) or fwhm_cur <= 0 or fwhm_cur > t_span
                    or not np.isfinite(t0_cur) or t0_cur < t_min - t_span or t0_cur > t_max
                    or (st['fw_wall'] and fwhm_cur < fw_lo)):
                return loss, A, fit_m
            if ca_free and (ca_w_cur < fwhm_cur / 5 or ca_w_cur > 5 * fwhm_cur):
                return loss, A, fit_m
        try:
            Kc = target_rate_matrix(sch, tau_cur)
            C_cur = target_conc(t, Kc, sch.c0, t0_cur, fwhm_cur, irf_mode)[0]
        except Exception:
            return loss, A, fit_m
        st['C'] = C_cur
        Cca = st['Cca']
        if ca_k > 0:
            Cca = ca_basis(t, t0_cur, ca_w_cur if ca_own else fwhm_cur, ca_k)
            st['Cca'] = Cca
        n_lin = 0
        for gg in range(n_g):
            Tm2 = Tg[gg]
            n_u = Tm2.shape[1]
            rws = np.flatnonzero(grp == gg)
            if n_u == 0 and ca_k == 0:
                row_done[rws] = True
                continue
            Bg = np.column_stack([C_cur @ Tm2, Cca])
            nnc = np.concatenate([np.asarray(nnG[gg], dtype=bool).ravel(),
                                  np.zeros(ca_k, dtype=bool)])
            fin = row_finite[rws]
            if fin.any():
                r1 = rws[fin]
                Bs = Bg[t_sel, :]
                Xa = lstsq_min_norm(Bs, D[np.ix_(r1, t_sel)].T).T
                if nnc.any():
                    tol_n = 1e-10 * np.max(np.abs(Xa))
                    bad = np.flatnonzero(np.any(Xa[:, nnc] < -tol_n, axis=1))
                    for ib in bad:
                        Xa[ib, :] = _nn_solve(Bs, D[r1[ib], t_sel], nnc)
                    Xa[:, nnc] = np.maximum(Xa[:, nnc], 0)
                A[r1, :] = Xa[:, :n_u] @ Tm2.T
                Aca[r1, :] = Xa[:, n_u:]
                row_done[r1] = True
                n_lin += int(np.count_nonzero(Xa))
            for r2 in rws[~fin]:
                mm = t_sel & np.isfinite(D[r2, :])
                if np.count_nonzero(mm) <= n_u + ca_k:
                    continue
                xa = lstsq_min_norm(Bg[mm, :], D[r2, mm])
                if nnc.any() and np.any(xa[nnc] < -1e-10 * np.max(np.abs(xa))):
                    xa = _nn_solve(Bg[mm, :], D[r2, mm], nnc)
                if nnc.any():
                    xa[nnc] = np.maximum(xa[nnc], 0)
                A[r2, :] = xa[:n_u] @ Tm2.T
                Aca[r2, :] = xa[n_u:]
                row_done[r2] = True
                n_lin += int(np.count_nonzero(xa))
        st['n_lin'] = n_lin
        if np.any(~np.isfinite(A)) or np.any(~np.isfinite(Aca)):
            A[:] = 0
            Aca[:] = 0
            return loss, A, fit_m
        fit_m = A @ C_cur.T + Aca @ Cca.T
        R = D - fit_m
        mfin = np.isfinite(R) & t_sel[None, :]
        st['n_loss'] = int(np.count_nonzero(mfin))
        if st['n_loss'] > 0:
            loss = float(np.sum(R[mfin] ** 2))
        return loss, A, fit_m

    # evaluate the model once outside the objective, so that a scheme that
    # cannot be computed reports why
    target_conc(t, target_rate_matrix(sch, tau_init), sch.c0, t0_init, fwhm_init, irf_mode)

    reordered = False
    if x0.size == 0:
        loss, A, fit_mat = objective(x0)
        iters, converged = 0, True
        init_loss = loss
        n_init = n_opt = st['n_loss']
        x_opt = x0
    else:
        init_loss = objective(x0)[0]
        n_init = st['n_loss']
        if init_loss >= 1e30:
            raise ValueError('The starting values are outside the allowed range (time constants between '
                             '1e-7 and 1e3 times the fitted delay span; t0 inside the data; FWHM below '
                             'the span).')
        kw = dict(tol_x=1e-8, tol_fun=1e-10, max_iter=5000, max_fun_evals=20000)
        x_opt, loss, flag, iters = fminsearch(lambda x: objective(x)[0], x0, **kw)
        converged = flag > 0
        # equivalent solutions: keep the ORDER of the starting values
        swp = np.flatnonzero(~has_lim)
        if swp.size > 1:
            rk_init = np.argsort(tau_init[idx_free[swp]], kind='stable')
            t_fit = np.exp(x_opt[swp])
            t_alt = t_fit.copy()
            t_alt[rk_init] = np.sort(t_fit)
            if np.any(np.abs(t_alt - t_fit) > 1e-9 * t_fit):
                x_alt = x_opt.copy()
                x_alt[swp] = np.log(t_alt)
                if objective(x_alt)[0] <= loss * (1 + 1e-3):
                    x_alt, loss_alt, _, _ = fminsearch(lambda x: objective(x)[0], x_alt, **kw)
                    if loss_alt <= loss * (1 + 1e-6):
                        x_opt, loss = x_alt, loss_alt
                        reordered = True
        _, A, fit_mat = objective(x_opt)
        n_opt = st['n_loss']
    n_lin_opt = st['n_lin']
    tau_out = st['tau'].copy()
    t0_out, fwhm_out = float(st['t0']), float(st['fwhm'])
    C_opt = st['C'].copy()

    at_limit = np.zeros(n_p, dtype=bool)
    for q in range(n_free):
        if has_lim[q]:
            fr0 = (np.log(tau_out[idx_free[q]]) - l_lo[q]) / (l_hi[q] - l_lo[q])
            at_limit[idx_free[q]] = fr0 < 1e-4 or fr0 > 1 - 1e-4

    st['fw_wall'] = False        # the search is over
    tau_std = np.full(n_p, np.nan)
    t0_std = np.nan
    fwhm_std = np.nan
    ca_std = np.nan
    if x0.size:
        use_u = ~at_limit[idx_free]
        p_opt = list(tau_out[idx_free[use_u]])
        n_tau = len(p_opt)
        if not t0_fixed:
            p_opt.append(t0_out)
        if not fwhm_fixed:
            p_opt.append(fwhm_out)
        if ca_free:
            p_opt.append(float(st['ca_w']))
        p_opt = np.asarray(p_opt, dtype=float)
        res_mask = np.isfinite(D) & t_sel[None, :]

        def residual_at_phys(pp):
            tv = tau_out[idx_free].copy()
            tv[use_u] = pp[:n_tau]
            xx = np.zeros(x0.size)
            xx[:n_free] = tau2x(tv)
            p3, p4 = n_tau, n_free
            if not t0_fixed:
                xx[p4] = pp[p3]
                p3 += 1
                p4 += 1
            if not fwhm_fixed:
                xx[p4] = np.log(max(pp[p3], _REALMIN))
                p3 += 1
                p4 += 1
            if ca_free:
                xx[p4] = np.log(max(pp[p3], _REALMIN))
            fm = objective(xx)[2]
            r = (D - fm)[res_mask]
            r[~np.isfinite(r)] = 0
            return r

        if p_opt.size:
            sig = fit_param_uncertainty(residual_at_phys, p_opt, loss, n_opt - n_lin_opt)
            tau_std[idx_free[use_u]] = sig[:n_tau]
            p2 = n_tau
            if not t0_fixed:
                t0_std = sig[p2]
                p2 += 1
            if not fwhm_fixed:
                fwhm_std = sig[p2]
                p2 += 1
            if ca_free:
                ca_std = sig[p2]
        _, A, fit_mat = objective(x_opt)

    sas = A.copy()
    sas[~st['row_done'], :] = np.nan
    ca_spec = st['Aca'].copy()
    ca_spec[~st['row_done'], :] = np.nan
    ca_w_out = np.nan
    if ca_k > 0:
        ca_w_out = float(st['ca_w']) if ca_own else fwhm_out
    K_opt = target_rate_matrix(sch, tau_out)
    _, obs_tau, mode_w = target_conc(t, K_opt, sch.c0, t0_out, fwhm_out, irf_mode)
    fin_tau = obs_tau[np.isfinite(obs_tau)]
    n_obs = 0
    if fin_tau.size:
        n_obs = 1 + int(np.count_nonzero(np.diff(fin_tau) > 1e-6 * fin_tau[:-1]))
    res = Result(tau=tau_out, t0=t0_out, fwhm=fwhm_out, sas=sas, fit=fit_mat,
                 rss=loss, iters=iters, converged=bool(converged),
                 rms=float(np.sqrt(loss / max(1, n_opt))),
                 initial_rms=float(np.sqrt(init_loss / max(1, n_init))),
                 n_used=n_opt, tau_std=tau_std, t0_std=t0_std, fwhm_std=fwhm_std,
                 at_limit=at_limit, conc=C_opt, K=K_opt, obs_tau=obs_tau, n_obs=n_obs,
                 das=sas @ mode_w, n_free=n_free, t_sel=t_sel, reordered=reordered,
                 scheme=sch, ca_order=ca_k, ca_spec=ca_spec, ca_fwhm=ca_w_out,
                 ca_fwhm_std=ca_std, ca_mode=ca['mode'],
                 fwhm_lo=np.nan, fwhm_at_limit=False, irf_retry=0)
    if not fwhm_fixed and fw_res > 0:
        res.fwhm_lo = fw_res
        res.fwhm_at_limit = irf_unresolved(fwhm_out, fw_res)
    if (not fwhm_fixed and fw_mode == 'auto' and x0.size
            and irf_width_lost(fwhm_out, fwhm_std, fw_res)):
        # the direct fit lost the IRF width: IRF fixed first, then all free
        fw_s0 = float(fwhm_init)
        if irf_unresolved(fwhm_init, fw_res):
            fw_s0 = 4 * fw_res
        try:
            s1 = fit_target_analysis(D, wl, t, sch, tau_init, tau_fixed, t0_init, fw_s0,
                                     True, True, irf_mode, t_range, tau_lo, tau_hi, ca, 'plain')
            tau_s1 = s1.tau.copy()
            q_lo = ~tau_fixed & np.isfinite(tau_lo) & (tau_lo > 0)
            q_hi = ~tau_fixed & np.isfinite(tau_hi) & (tau_hi > 0)
            tau_s1[q_lo] = np.maximum(tau_s1[q_lo], tau_lo[q_lo])
            tau_s1[q_hi] = np.minimum(tau_s1[q_hi], tau_hi[q_hi])
            s2 = fit_target_analysis(D, wl, t, sch, tau_s1, tau_fixed, t0_init, fw_s0,
                                     t0_fixed_arg, fwhm_fixed_arg, irf_mode, t_range,
                                     tau_lo, tau_hi, ca, 'limit')
            it_all = res.iters + s1.iters + s2.iters
            if (s2.rss < res.rss * (1 - 1e-6)
                    or (res.fwhm_at_limit and s2.rss <= res.rss * (1 + 1e-6))):
                s2.initial_rms = res.initial_rms
                res = s2
                res.irf_retry = 1
            else:
                res.irf_retry = 2
            res.iters = it_all
        except Exception:
            pass
    return res


def target_step_table(sch, tau):
    """Per step: label, rate expression, time constant of the step
    (tau of its parameter / multiplier) and the branching fraction
    k / (sum of k leaving the same compartment)."""
    tau = np.asarray(tau, dtype=float).ravel()
    k = np.array([s['mult'] / tau[s['par']] for s in sch.steps])
    rows = []
    for i, s in enumerate(sch.steps):
        tot = sum(k[j] for j, q in enumerate(sch.steps) if q['from'] == s['from'])
        rows.append({'label': s['label'], 'expr': s['expr'], 'tau': 1.0 / k[i],
                     'branching': k[i] / tot if tot > 0 else np.nan})
    return rows
