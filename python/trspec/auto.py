"""Automatic analysis: pre-processing (noise, time zero, background, mask,
chirp, crop), model choice and report.

Port of autoOptions ... autoReportLines of TRSpecAnalyzer.m. Positions
along the delay or wavelength axis that the MATLAB code keeps as indices
(onset 'on', 'jOn', window limits) are kept 1-based here as well, so that
the two codes can be read side by side; arrays are indexed with value - 1.
"""
import math

import numpy as np

from ._util import Result, interp1_nan
from .chirp import chirp_model, fit_chirp_params
from .kinetics import matlab_round

_REALMIN = np.finfo(float).tiny


def auto_options(o=None):
    """Options of the automatic analysis with their defaults."""
    d = dict(doBg=True, doMask=True, doChirp=True, doCrop=True, maxComp=5,
             minGain=0.05, zOnset=4, maskK=3, maskSN=10, maskSevere=10,
             nChirpPts=10, ca='auto', caGain=0.15, caTauMin=1.0 / 3,
             tryAll=False)
    out = dict(o) if isinstance(o, dict) else {}
    for k, v in d.items():
        if k not in out or out[k] is None or \
                (not isinstance(out[k], (str, bool)) and np.size(out[k]) == 0):
            out[k] = v
    out['maxComp'] = int(max(1, min(5, matlab_round(out['maxComp']))))
    if not isinstance(out['ca'], str) or out['ca'] not in ('auto', 'off', 'on'):
        out['ca'] = 'auto'
    out['tryAll'] = bool(np.ravel(out['tryAll'])[0])
    return out


def _nanmean(X, axis):
    X = np.asarray(X, float)
    ok = np.isfinite(X)
    n = ok.sum(axis=axis)
    with np.errstate(invalid='ignore', divide='ignore'):
        m = np.where(ok, X, 0.0).sum(axis=axis) / n
    return np.where(n == 0, np.nan, m)


def _median(x):
    x = np.asarray(x, float)
    return float(np.median(x)) if x.size else np.nan


def _run_med(x, w):
    """Running median over w points, NaN ignored."""
    x = np.asarray(x, float).ravel()
    n = x.size
    y = np.full(n, np.nan)
    h = int(math.floor(w / 2))
    for i in range(n):
        seg = x[max(0, i - h):min(n, i + h + 1)]
        seg = seg[np.isfinite(seg)]
        if seg.size:
            y[i] = np.median(seg)
    return y


def auto_ac1(x):
    """Lag-1 autocorrelation of the finite values of x."""
    x = np.asarray(x, float).ravel()
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    x = x - np.mean(x)
    den = np.sum(x ** 2)
    if den > 0:
        return float(np.sum(x[:-1] * x[1:]) / den)
    return np.nan


def auto_noise(D, t):
    """Noise of every wavelength from the deviation of each point from the
    line through its two neighbors (robust)."""
    D = np.asarray(D, float)
    t = np.asarray(t, float).ravel()
    M = D.shape[0]
    sg = np.full(M, np.nan)
    for i in range(M):
        m = np.isfinite(D[i, :])
        if np.count_nonzero(m) < 8:
            continue
        y = D[i, m]
        x = t[m]
        with np.errstate(invalid='ignore', divide='ignore'):
            w1 = (x[2:] - x[1:-1]) / (x[2:] - x[:-2])
            w2 = 1 - w1
            r = y[1:-1] - (w1 * y[:-2] + w2 * y[2:])
            r = r / np.sqrt(1 + w1 ** 2 + w2 ** 2)
        r = r[np.isfinite(r)]
        if r.size < 5:
            continue
        sg[i] = 1.4826 * np.median(np.abs(r - np.median(r)))
    return sg


def auto_rise(D, t, sg, n_base, z_thr):
    """Where the signal of every wavelength leaves its baseline.
    Returns a Result with b (baseline), on, r30, r80 (1-based delay
    indices, NaN if none) and abrupt."""
    D = np.asarray(D, float)
    M, N = D.shape
    t = np.asarray(t, float).ravel()
    n_base = int(n_base)
    b = np.full(M, np.nan)
    for i in range(M):
        seg = D[i, :n_base]
        seg = seg[np.isfinite(seg)]
        if seg.size:
            b[i] = np.mean(seg)
    thr = z_thr * sg * np.sqrt(1 + 1.0 / n_base)
    on = np.full(M, np.nan)
    r30 = np.full(M, np.nan)
    r80 = np.full(M, np.nan)
    abrupt = np.zeros(M, bool)
    for i in range(M):
        if not np.isfinite(sg[i]) or sg[i] <= 0 or not np.isfinite(b[i]):
            continue
        d = D[i, :] - b[i]
        with np.errstate(invalid='ignore'):
            big = np.abs(d) > thr[i]
        j = n_base + 1                      # 1-based
        found = False
        while j < N:
            if big[j - 1] and big[j] and d[j - 1] * d[j] > 0:
                found = True
                break
            j += 1
        if not found:
            continue
        on[i] = j
        ad = np.abs(d)
        ad[~np.isfinite(ad)] = 0
        tail = ad[j - 1:]
        mx = np.max(tail)
        k = int(np.argmax(tail >= 0.3 * mx)) + j          # 1-based
        r30[i] = k
        r80[i] = int(np.argmax(tail >= 0.8 * mx)) + j
        abrupt[i] = (k - j) <= 4 or \
            (t[k - 1] - t[j - 1]) <= 0.02 * (t[-1] - t[j - 1])
    return Result(b=b, on=on, r30=r30, r80=r80, abrupt=abrupt)


def auto_ridge(D, t, sg, jA, jB, w):
    """Time of the steepest change of every wavelength between the delays
    jA and jB (1-based), from differences over w points.
    Returns (tr, zj, wid, dtk)."""
    D = np.asarray(D, float)
    t = np.asarray(t, float).ravel()
    M = D.shape[0]
    tr = np.full(M, np.nan)
    zj = np.zeros(M)
    wid = np.full(M, np.nan)
    dtk = np.full(M, np.nan)
    jA, jB, w = int(jA), int(jB), int(w)
    ts = t[jA - 1:jB]
    nR = ts.size - w
    if nR < 3:
        return tr, zj, wid, dtk
    dt = ts[w:] - ts[:-w]
    tm = (ts[w:] + ts[:-w]) / 2
    d1 = np.diff(ts)
    for i in range(M):
        if not np.isfinite(sg[i]) or sg[i] <= 0:
            continue
        y = D[i, jA - 1:jB]
        dA = y[w:] - y[:-w]
        with np.errstate(invalid='ignore', divide='ignore'):
            rate = np.abs(dA / dt)
        rate[~np.isfinite(rate)] = 0
        k = int(np.argmax(rate))                 # 0-based from here on
        mx = rate[k]
        if not (mx > 0):
            continue
        zj[i] = abs(dA[k]) / (sg[i] * np.sqrt(2))
        a = k
        while a > 0 and rate[a - 1] >= mx / 2:
            a -= 1
        b = k
        while b < nR - 1 and rate[b + 1] >= mx / 2:
            b += 1
        dtk[i] = d1[min(k, d1.size - 1)]
        tr[i] = tm[k]
        wLobe = tm[b] - tm[a] + dtk[i]
        a3 = k
        while a3 > 0 and rate[a3 - 1] >= 0.3 * mx:
            a3 -= 1
        b3 = k
        while b3 < nR - 1 and rate[b3 + 1] >= 0.3 * mx:
            b3 += 1
        with np.errstate(invalid='ignore'):
            sel = (np.abs(tm - tm[k]) <= 3 * wLobe) & (rate >= 0.3 * mx) & \
                (np.abs(dA) >= 4 * np.sqrt(2) * sg[i])
        sel[a3:b3 + 1] = True
        several = bool(np.any(sel[:a3]) or np.any(sel[b3 + 1:]))
        if several:
            wgt = rate[sel] ** 2
            tr[i] = np.sum(wgt * tm[sel]) / np.sum(wgt)
        elif b - a + 1 > 4:
            tr[i] = np.sum(rate[a:b + 1] * tm[a:b + 1]) / np.sum(rate[a:b + 1])
        elif 0 < k < nR - 1 and \
                abs((tm[k] - tm[k - 1]) - (tm[k + 1] - tm[k])) < \
                0.01 * (tm[k + 1] - tm[k]):
            den = rate[k - 1] - 2 * rate[k] + rate[k + 1]
            if den < 0:
                tr[i] = tm[k] + 0.5 * (rate[k - 1] - rate[k + 1]) / den * \
                    (tm[k + 1] - tm[k])
        wRise = np.sqrt(max(wLobe ** 2 - 0.46 * dt[k] ** 2, dtk[i] ** 2))
        iS = np.flatnonzero(sel)
        ext = tm[iS[-1]] - tm[iS[0]] + dtk[i]
        wid[i] = max(wRise, ext / 1.9)
    return tr, zj, wid, dtk


def _runs(flag):
    """Starts and ends (1-based, inclusive) of the runs of True."""
    f = np.asarray(flag, bool).astype(int)
    d = np.diff(np.concatenate([[0], f, [0]]))
    return np.flatnonzero(d == 1) + 1, np.flatnonzero(d == -1)


def _mad(x):
    x = np.asarray(x, float)
    if x.size == 0:
        return np.nan
    return 1.4826 * float(np.median(np.abs(x - np.median(x))))


def auto_preprocess(wl, t, D, opt=None, time_unit='ps', spec_unit='nm'):
    """What the data need before a fit: noise of every wavelength, time
    zero, background, noisy channels, chirp, crop. D is n_wavelength x
    n_delay; time_unit / spec_unit are the unit names used in messages.

    Returns a Result P (fields as in the MATLAB program: opt, noise,
    rowOK, hasT0, nBg, bg, bgStatus, maskRows, maskReg, keptReg, chirp,
    keepT, keepW, cropT, cropW, Dfull, wl, t, D, t0Guess, fwGuess, changed,
    checks, ...)."""
    opt = auto_options(opt)
    wl = np.asarray(wl, float).ravel()
    t = np.asarray(t, float).ravel()
    D = np.asarray(D, float)
    M, N = D.shape
    if M < 5 or N < 10:
        raise ValueError('The data set is too small for the automatic '
                         'analysis (%d x %d).' % (M, N))
    tu, xu = time_unit, spec_unit
    checks = []
    sg = auto_noise(D, t)
    with np.errstate(invalid='ignore'):
        rowOK = np.isfinite(sg) & (sg > 0)
    if np.count_nonzero(rowOK) < 5:
        raise ValueError('Fewer than 5 wavelengths hold usable data.')

    # ---- onset of the signal, number of baseline delays --------------
    nBase, jOn, nBg = 3, np.nan, 0
    for _ in range(6):
        R = auto_rise(D, t, sg, nBase, opt['zOnset'])
        det = rowOK & np.isfinite(R.on)
        if np.count_nonzero(det) < max(3, math.ceil(0.05 * np.count_nonzero(rowOK))):
            jOn, nBg = np.nan, 0
            break
        onS = np.sort(R.on[det])
        jOn = int(onS[max(1, math.ceil(0.05 * onS.size)) - 1])
        nBg = jOn - 1 - max(2, math.ceil(0.15 * (jOn - 1)))
        if nBg <= nBase:
            break
        nBase = nBg
    R = auto_rise(D, t, sg, nBase, opt['zOnset'])
    det = rowOK & np.isfinite(R.on)
    hasT0 = bool(np.isfinite(jOn))
    t0Why = ''
    if not hasT0:
        t0Why = 'no signal onset was found'
    elif np.count_nonzero(det & R.abrupt) < 0.2 * np.count_nonzero(det):
        hasT0 = False
        t0Why = 'the signal changes only gradually after the first delays'

    # ---- signal size of every wavelength ------------------------------
    b0 = np.where(np.isfinite(R.b), R.b, 0.0)
    Xa = np.abs(D - b0[:, None])
    Xa[~np.isfinite(Xa)] = 0
    if N >= 3:
        X3 = np.sort(np.stack([Xa[:, :N - 2], Xa[:, 1:N - 1], Xa[:, 2:N]],
                              axis=2), axis=2)
        sigAmp = np.max(X3[:, :, 1], axis=1)
    else:
        sigAmp = np.max(Xa, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        snr = sigAmp / sg
    snr[~rowOK] = np.nan

    # ---- noisy channels -------------------------------------------------
    maskRows = np.zeros(M, bool)
    maskStatus = 'off'
    ratio = np.full(M, np.nan)
    noisy = np.zeros(M, bool)
    if opt['doMask']:
        maskStatus = 'none'
        w = max(11, 2 * int(matlab_round(0.1 * M)) + 1)
        sgf = np.where(rowOK, sg, np.nan)
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = sg / _run_med(sgf, w)
            noisy = rowOK & (ratio > opt['maskK'])
            sev = (snr < opt['maskSN']) | (ratio > opt['maskSevere'])
            sevW = (snr < opt['maskSN']) | (ratio > opt['maskSevere'] / 2)
        flag = noisy & sev
        if np.any(flag):
            lowK = max(1.5, opt['maskK'] / 2)
            fl = flag.copy()
            with np.errstate(invalid='ignore'):
                grow = rowOK & (ratio > lowK) & sevW
            for q in range(1, M):
                if fl[q - 1] and grow[q]:
                    fl[q] = True
            for q in range(M - 2, -1, -1):
                if fl[q + 1] and grow[q]:
                    fl[q] = True
            flag = fl.copy()
            fl[:-1] = fl[:-1] | flag[1:]
            fl[1:] = fl[1:] | flag[:-1]
            idx = np.flatnonzero(fl)
            for q in range(idx.size - 1):
                if idx[q + 1] - idx[q] <= 3:
                    fl[idx[q]:idx[q + 1] + 1] = True
            if np.count_nonzero(fl) > 0.25 * M:
                maskStatus = 'toomany'
            else:
                maskStatus = 'applied'
                maskRows = fl
    rs, re = _runs(maskRows)
    edgeLo = edgeHi = 0
    if opt['doCrop'] and rs.size:
        if rs[0] == 1:
            edgeLo = int(re[0])
        if re[-1] == M and (rs.size > 1 or edgeLo == 0):
            edgeHi = int(rs[-1])

    def region(a, b):
        seg = slice(a - 1, b)
        return [float(np.min(wl[seg])), float(np.max(wl[seg])),
                float(np.max(ratio[seg])), float(np.min(snr[seg]))]

    maskReg = []
    for q in range(rs.size):
        if (edgeLo > 0 and rs[q] == 1) or (edgeHi > 0 and re[q] == M):
            continue
        maskReg.append(region(int(rs[q]), int(re[q])))
    kn = noisy & ~maskRows
    ks, ke = _runs(kn)
    keptReg = [region(int(ks[q]), int(ke[q])) for q in range(ks.size)]
    maskReg = np.asarray(maskReg, float).reshape(-1, 4)
    keptReg = np.asarray(keptReg, float).reshape(-1, 4)
    vr = rowOK & ~maskRows

    # ---- time zero along the spectrum (steepest rise) -------------------
    ch = dict(status='off', msg='', pts=np.zeros((0, 2)), p=None, rms=np.nan,
              t0=np.zeros(M), span=np.nan, fracIn=np.nan, scatter=np.nan,
              tol=np.nan, ridgeT=np.full(M, np.nan), good=np.zeros(M, bool),
              win=[np.nan, np.nan], dtZ=np.nan, wStep=1,
              cover=[np.nan, np.nan], model='', shift=0.0, monotonic=True,
              nRejected=0)
    fwGuess = t0Guess = np.nan
    tr = np.full(M, np.nan)
    good = np.zeros(M, bool)
    if hasT0:
        use = vr & np.isfinite(R.on) & R.abrupt
        if np.count_nonzero(use) >= 5:
            onU = np.sort(R.on[use])
            j0 = int(onU[max(1, int(matlab_round(0.01 * onU.size))) - 1])
            r8 = np.sort(R.r80[use])
            j1 = int(r8[max(1, math.ceil(0.98 * r8.size)) - 1])
            wList = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]
            nNeed = max(5, math.ceil(0.3 * np.count_nonzero(use)))
            bestN = -1
            zj = wid = dtk = None
            for q, ws in enumerate(wList):
                jA = max(1, j0 - 5 - ws)
                jB = min(N, j1 + 3 + ws)
                if jB - jA < 4 * ws:
                    if q > 0:
                        break
                    jB = min(N, jA + 4)
                trQ, zjQ, widQ, dtkQ = auto_ridge(D, t, sg, jA, jB, ws)
                nClear = int(np.count_nonzero(use & (zjQ > 6)))
                if nClear > bestN:
                    bestN = nClear
                    tr, zj, wid, dtk = trQ, zjQ, widQ, dtkQ
                    ch['win'] = [float(t[jA - 1]), float(t[jB - 1])]
                    ch['wStep'] = ws
                if nClear >= nNeed:
                    break
            good = use & np.isfinite(tr) & (zj > 5)
            ch['ridgeT'] = tr
            ch['good'] = good
            if np.count_nonzero(good) >= 5:
                ch['dtZ'] = _median(dtk[good])
                fwGuess = _median(wid[good])
                t0Guess = _median(tr[good])
        if np.count_nonzero(good) < max(5, math.ceil(0.05 * np.count_nonzero(vr))) \
                or (nBg < 3 and ch['wStep'] > 1):
            hasT0 = False
            t0Why = 'no sharp rise: the signal seems to be there from the ' \
                    'first delays on'
            fwGuess = t0Guess = np.nan
            good = np.zeros(M, bool)
            ch['good'] = good

    # ---- background ------------------------------------------------------
    bg = np.zeros(M)
    bgStatus = 'off'
    bSig = np.nan
    if opt['doBg']:
        if not hasT0:
            bgStatus = 'notime0'
        elif nBg < 3:
            bgStatus = 'few'
        else:
            b = _nanmean(D[:, :nBg], 1)
            e = _nanmean(D[:, max(1, N - 2) - 1:N], 1)
            m = rowOK & np.isfinite(b) & np.isfinite(e)
            bw = b[m] / sg[m]
            ew = e[m] / sg[m]
            dw = ew - bw
            nb = float(np.linalg.norm(bw))
            bSig = nb / np.sqrt(max(1, np.count_nonzero(m))) * np.sqrt(nBg)
            cosv = float(bw @ dw) / max(nb * float(np.linalg.norm(dw)), _REALMIN)
            if bSig > 3 and np.linalg.norm(ew) < 0.5 * nb and cosv < -0.8:
                bgStatus = 'signal'
            elif bSig < 2:
                bgStatus = 'zero'
            else:
                bgStatus = 'applied'
                bg = np.where(np.isfinite(b), b, 0.0)
    Db = D - bg[:, None]

    # ---- chirp -------------------------------------------------------------
    if not hasT0:
        ch['status'] = 'notime0'
    else:
        nB = min(int(opt['nChirpPts']), np.count_nonzero(good) // 3)
        pts = []
        if nB >= 5:
            edges = np.linspace(np.min(wl[vr]), np.max(wl[vr]), nB + 1)
            for k in range(nB):
                if k < nB - 1:
                    inb = (wl >= edges[k]) & (wl < edges[k + 1])
                else:
                    inb = (wl >= edges[k]) & (wl <= edges[k + 1])
                g = inb & good
                if np.count_nonzero(g) >= max(3, math.ceil(0.3 * np.count_nonzero(inb & vr))):
                    pts.append([np.median(wl[g]), np.median(tr[g])])
            pts = np.asarray(pts, float).reshape(-1, 2)
        else:
            pts = np.column_stack([wl[good], tr[good]])
        ch['tol'] = float(np.nanmax([2 * ch['dtZ'], 0.3 * fwGuess])) \
            if np.isfinite(ch['dtZ']) or np.isfinite(fwGuess) else np.nan
        if pts.shape[0] < 5:
            ch['status'] = 'failed'
            ch['msg'] = 'fewer than 5 usable points along the spectrum'
        else:
            for attempt in (1, 2):
                try:
                    if attempt == 1:
                        p, rmsF = fit_chirp_params(pts)
                        r = tr - chirp_model(p, wl)
                        rg = r[good]
                        s = _mad(rg)
                        with np.errstate(invalid='ignore'):
                            inl = good & (np.abs(r - np.median(rg))
                                          <= max(3 * s, 0.75 * ch['tol']))
                        if np.count_nonzero(inl) >= 8:
                            pts = np.column_stack([wl[inl], tr[inl]])
                            p, rmsF = fit_chirp_params(pts)
                        t0m = chirp_model(p, wl)
                        ch['model'] = 'formula'
                    else:
                        xm = float(np.mean(wl[good]))
                        xs = float(np.std(wl[good], ddof=1))
                        xg = (wl[good] - xm) / xs
                        Ap = np.column_stack([np.ones(xg.size), xg, xg ** 2])
                        cp = np.linalg.lstsq(Ap, tr[good], rcond=None)[0]
                        rg = tr[good] - Ap @ cp
                        s = _mad(rg)
                        keepP = np.abs(rg - np.median(rg)) <= \
                            max(3 * s, 0.75 * ch['tol'])
                        if np.count_nonzero(keepP) < 8:
                            keepP = np.ones(rg.size, bool)
                        gi = np.flatnonzero(good)
                        pts = np.column_stack([wl[gi[keepP]], tr[gi[keepP]]])
                        cp = np.linalg.lstsq(Ap[keepP, :], pts[:, 1],
                                             rcond=None)[0]
                        xa = (wl - xm) / xs
                        t0m = np.column_stack([np.ones(M), xa, xa ** 2]) @ cp
                        p = None
                        o = np.argsort(wl)
                        rmsF = float(np.sqrt(np.mean(
                            (pts[:, 1] - interp1_nan(wl[o], t0m[o],
                                                     pts[:, 0])) ** 2)))
                        ch['model'] = 'parabola'
                    if np.any(~np.isfinite(t0m[vr])):
                        raise ValueError('the fitted curve is not defined '
                                         'over the whole spectrum')
                    lo, hi = float(np.min(pts[:, 1])), float(np.max(pts[:, 1]))
                    slack = 0.5 * (hi - lo) + ch['tol']
                    if np.any((t0m[vr] < lo - slack) | (t0m[vr] > hi + slack)):
                        raise ValueError(
                            'the fitted curve runs far outside the measured '
                            'rises (%.3g ... %.3g %s)' % (lo, hi, tu))
                    res = tr[good] - t0m[good]
                    ch['nRejected'] = int(np.count_nonzero(good) - pts.shape[0])
                    ch['cover'] = [float(np.min(pts[:, 0])),
                                   float(np.max(pts[:, 0]))]
                    ch['pts'] = pts
                    ch['p'] = p
                    ch['rms'] = float(rmsF)
                    ch['fracIn'] = float(np.count_nonzero(
                        np.abs(res) <= ch['tol']) / res.size)
                    ch['scatter'] = _mad(res)
                    ch['span'] = float(np.max(t0m[vr]) - np.min(t0m[vr]))
                    d1 = np.diff(t0m[vr])
                    ch['monotonic'] = bool(np.all(d1 >= 0) or np.all(d1 <= 0))
                    inC = vr & (wl >= ch['cover'][0]) & (wl <= ch['cover'][1])
                    spanC = float(np.max(t0m[inC]) - np.min(t0m[inC]))
                    if spanC < 0.75 * ch['tol']:
                        ch['status'] = 'none'
                        ch['span'] = spanC
                    elif ch['fracIn'] < 0.7:
                        raise ValueError(
                            'only %.0f %% of the wavelengths lie within '
                            '%.2g %s of the fitted curve'
                            % (100 * ch['fracIn'], ch['tol'], tu))
                    elif not opt['doChirp']:
                        ch['status'] = 'found'
                        ch['t0'] = t0m
                    else:
                        ch['status'] = 'applied'
                        ch['t0'] = t0m
                    ch['msg'] = ''
                    break
                except Exception as errC:
                    ch['status'] = 'failed'
                    if attempt == 1:
                        ch['msg'] = str(errC)

    # ---- carry out: chirp, mask, crop ---------------------------------------
    Dc = Db.copy()
    tOut = t.copy()
    if ch['status'] == 'applied':
        try:
            for i in range(M):
                if not np.isfinite(ch['t0'][i]):
                    Dc[i, :] = np.nan
                    continue
                Dc[i, :] = interp1_nan(t, Db[i, :], t + ch['t0'][i])
            t0Guess = 0.0
        except Exception as errI:
            Dc = Db.copy()
            ch['status'] = 'failed'
            ch['msg'] = 'the correction could not be carried out: %s' % errI
    elif ch['status'] == 'none':
        if opt['doChirp'] and np.isfinite(t0Guess) and \
                abs(t0Guess) > 0.5 * ch['tol']:
            ch['shift'] = float(t0Guess)
            tOut = t - ch['shift']
            t0Guess = 0.0
    Dc[maskRows, :] = np.nan
    keepT = np.ones(N, bool)
    keepW = np.ones(M, bool)
    if opt['doCrop']:
        colOK = np.all(np.isfinite(Dc[vr, :]), axis=0)
        if np.any(colOK) and not np.all(colOK):
            c = np.flatnonzero(colOK)
            j1c, j2c = int(c[0]), int(c[-1])
            if j2c - j1c + 1 >= 10:
                keepT[:j1c] = False
                keepT[j2c + 1:] = False
        if edgeLo > 0:
            keepW[:edgeLo] = False
        if edgeHi > 0:
            keepW[edgeHi - 1:] = False

    # ---- things to check -----------------------------------------------------
    if not hasT0:
        checks.append('No time zero was found in the data (%s). Background '
                      'and chirp were left alone, and the fits use no '
                      'instrument response.' % t0Why)
    if bgStatus == 'few':
        checks.append('Only %d delay(s) lie safely before the signal: no '
                      'background was subtracted. Set it by hand if there is '
                      'one.' % max(nBg, 0))
    elif bgStatus == 'signal':
        checks.append('The first delays look like a signal that decays '
                      'later, not like a background: nothing was subtracted.')
    if maskStatus == 'toomany':
        checks.append('More than a quarter of the channels are much noisier '
                      'than their neighbors: no mask was set. Look at the '
                      'noise plot.')
    if keptReg.shape[0]:
        checks.append('Some channels are noisier than their neighbors but '
                      'were kept, because their signal is far above that '
                      'noise (listed under Mask). The fits are not weighted, '
                      'so such channels count more than quiet ones. If the '
                      'residual shows stripes there, mask them: Mask rule '
                      '"Every noisy channel", or Mask in the main window.')
    if ch['status'] == 'failed':
        checks.append('The time zero along the spectrum could not be '
                      'determined reliably (%s): no chirp correction. Use '
                      'the Chirp Corr window.' % ch['msg'])
    elif ch['status'] in ('applied', 'found'):
        if ch['model'] == 'parabola':
            checks.append('The chirp formula of the Chirp Corr window gave '
                          'no sensible curve here; a parabola through the '
                          'rises was used instead. Compare the curve with '
                          'the map.')
        if not ch['monotonic']:
            checks.append('The fitted time-zero curve is not monotonic along '
                          'the spectrum. Chirp from material dispersion is: '
                          'compare the curve with the map.')
        if ch['fracIn'] < 0.9:
            checks.append('Only %.0f %% of the wavelengths lie within %.2g '
                          '%s of the time-zero curve: compare the curve with '
                          'the map.' % (100 * ch['fracIn'], ch['tol'], tu))
        wr = float(np.max(wl[vr]) - np.min(wl[vr]))
        if ch['cover'][0] - np.min(wl[vr]) > 0.1 * wr or \
                np.max(wl[vr]) - ch['cover'][1] > 0.1 * wr:
            checks.append('A sharp rise is seen only from %.4g to %.4g %s. '
                          'Outside, the time-zero curve is an extrapolation.'
                          % (ch['cover'][0], ch['cover'][1], xu))

    P = Result()
    P.opt = opt
    P.wlIn, P.tIn = wl, t
    P.noise, P.rowOK, P.noiseRatio = sg, rowOK, ratio
    P.hasT0, P.t0Why = hasT0, t0Why
    P.jOn = jOn
    P.tOn = float(t[int(jOn) - 1]) if hasT0 else np.nan
    P.nBg, P.bg, P.bgStatus, P.bgSig = nBg, bg, bgStatus, bSig
    P.maskRows, P.maskReg, P.maskStatus = maskRows, maskReg, maskStatus
    P.keptReg, P.snr, P.noisy = keptReg, snr, noisy
    P.chirp = ch
    P.keepT, P.keepW = keepT, keepW
    P.cropT = P.cropW = None
    if np.any(~keepT):
        c = np.flatnonzero(keepT)
        P.cropT = [float(t[c[0]]), float(t[c[-1]])]
    if np.any(~keepW):
        P.cropW = [float(np.min(wl[keepW])), float(np.max(wl[keepW]))]
    P.Dfull = Dc
    P.wl, P.t, P.D = wl[keepW], tOut[keepT], Dc[np.ix_(keepW, keepT)]
    P.noiseOut = sg[keepW]
    P.t0Guess, P.fwGuess = t0Guess, fwGuess
    P.changed = bool(bgStatus == 'applied' or np.any(maskRows)
                     or ch['status'] == 'applied' or ch['shift'] != 0
                     or np.any(~keepT) or np.any(~keepW))
    P.checks = checks
    P.tu, P.xu = tu, xu
    return P


# ======================================================================
# model choice
# ======================================================================
def _logspace(a, b, n):
    """MATLAB logspace(a, b, n): for n = 1 the upper end."""
    n = int(n)
    if n <= 0:
        return np.zeros(0)
    if n == 1:
        return np.array([10.0 ** b])
    return np.logspace(a, b, n)


def auto_fit_ga(D, t, tau0, t0s, fws, has_inf, irf, fw_fixed, ca_k=0,
                irf_fixed=False):
    """One global fit of the automatic analysis (plain search). Returns a
    Result F with the time constants in ascending order, or None if the
    fit failed."""
    from .global_analysis import fit_global_analysis
    is_none = irf == 'none'
    try:
        o = fit_global_analysis(
            D, t, tau0, t0s, fws, False, bool(is_none or irf_fixed),
            bool(is_none or fw_fixed or irf_fixed), bool(has_inf), None, None,
            irf, None, {'order': int(ca_k), 'mode': 'irf'}, 'plain')
    except Exception:
        return None
    if not np.isfinite(o.rss) or o.rss >= 1e29:
        return None
    tau = np.atleast_1d(np.asarray(o.tau, float)).ravel()
    ix = np.argsort(tau, kind='stable')
    tau = tau[ix]
    nT = tau.size
    perm = np.concatenate([ix, [nT]]) if has_inf else ix
    A = np.asarray(o.dads, float)
    ca_order = int(o.ca_order)
    ca_spec = np.asarray(o.ca_spec, float).reshape(A.shape[0], -1) \
        if ca_order > 0 else np.zeros((A.shape[0], 0))
    return Result(tau=tau, tauStd=np.atleast_1d(np.asarray(o.tau_std, float))[ix],
                  t0=float(o.t0), fwhm=float(o.fwhm), t0Std=float(o.t0_std),
                  fwhmStd=float(o.fwhm_std), A=A[:, perm],
                  fit=np.asarray(o.fit, float), rss=float(o.rss),
                  rms=float(o.rms), hasInf=bool(has_inf), caK=ca_order,
                  caSpec=ca_spec, irfFixed=bool(irf_fixed and not is_none))


def _undet_long(F, tEnd):
    with np.errstate(invalid='ignore', divide='ignore'):
        relE = F.tauStd / F.tau
        undet = ~np.isfinite(relE) | (relE > 0.5)
    longT = (F.tau > 2 * tEnd) | ((F.tau > 0.5 * tEnd) & undet)
    return undet, longT


def auto_fit_sequence(Dv, Dk, tv, C, opt, ca_k=0, prog=None, p0=0.0, p1=1.0):
    """Global fits with 1, 2, ... components until one more component no
    longer pays. C: dict with irf, t0s, fws, dtZ, tauLo, tauHi, tEnd, tu
    (and 'seed' for the pass with artifact terms). Returns a Result T
    with the lists N, Rms, Gain, Inf, Tau, Std, Why, F and nStar, stopWhy,
    fwFixed, fwNote, caK."""
    irf, dtZ, tEnd = C['irf'], C['dtZ'], C['tEnd']
    tauLo, tauHi = C['tauLo'], C['tauHi']
    nMax = int(opt['maxComp'])
    tabN, tabRms, tabGain, tabInf = [], [], [], []
    tabTau, tabStd, tabWhy, tabF = [], [], [], []
    fwFixed, fwNote = False, ''
    useInf = False
    prevTau = None
    t0c, fwc = C['t0s'], C['fws']
    nStar, stopWhy, decided = 0, '', False
    seed = C.get('seed') if ca_k > 0 else None

    def put(n, F, gain, why):
        for lst, v in ((tabN, n), (tabRms, F.rms), (tabGain, gain),
                       (tabInf, bool(F.hasInf)), (tabTau, F.tau),
                       (tabStd, F.tauStd), (tabWhy, why), (tabF, F)):
            if len(lst) >= n:
                lst[n - 1] = v
            else:
                lst.append(v)

    n = 1
    while n <= nMax:
        if prog is not None:
            if ca_k > 0:
                msg = 'Global fit with %d component(s) and artifact terms...' % n
            else:
                msg = 'Global fit with %d component(s)...' % n
            prog(p0 + (p1 - p0) * (n - 1) / nMax, msg)
        below = seed is not None and n < seed['n']
        if below:
            F = auto_fit_ga(Dv, tv, seed['tau0'][n - 1], seed['F'].t0,
                            seed['F'].fwhm, bool(seed['inf0'][n - 1]), irf,
                            True, ca_k, True)
            if F is None:
                stopWhy = 'the fit with %d components failed' % n
                break
            gain = np.nan
            if n > 1:
                gain = 1 - F.rms / tabRms[n - 2]
            put(n, F, gain, '')
            n += 1
            continue
        atSeed = seed is not None and n == seed['n']
        if atSeed:
            useInf = bool(seed['F'].hasInf)
            fwFixed, fwNote = seed['fwFixed'], seed['fwNote']
        m = n - int(useInf)
        st = []
        if atSeed:
            t0c, fwc = seed['F'].t0, seed['F'].fwhm
        elif prevTau is None or np.size(prevTau) == 0:
            st.append(np.array([max(tEnd / 5, 2 * tauLo)]))
            st.append(np.array([np.sqrt(2 * tauLo * tauHi)]))
        else:
            st.append(_logspace(np.log10(2 * tauLo), np.log10(tauHi), m))
            pv = np.sort(np.asarray(prevTau, float))
            gr = np.concatenate([[tauLo], pv, [max(2 * tEnd, 2 * pv[-1])]])
            g = int(np.argmax(np.diff(np.log(gr))))
            st.append(np.sort(np.concatenate([pv, [np.sqrt(gr[g] * gr[g + 1])]])))
            if pv[0] / 4 > tauLo:
                st.append(np.concatenate([[pv[0] / 4], pv]))
            if pv[-1] * 4 < 3 * tEnd:
                st.append(np.concatenate([pv, [pv[-1] * 4]]))
        best = None
        for s in st:
            Fs = auto_fit_ga(Dk, tv, s, t0c, fwc, useInf, irf, fwFixed, ca_k)
            if Fs is not None and (best is None or Fs.rss < best.rss):
                best = Fs
        F = None
        if atSeed:
            F = seed['F']
        elif best is not None:
            F = auto_fit_ga(Dv, tv, best.tau, best.t0, best.fwhm, useInf, irf,
                            fwFixed, ca_k)
        if F is None:
            if not decided:
                stopWhy = 'the fit with %d components failed' % n
            break
        if irf == 'numerical' and not fwFixed and not atSeed and \
                F.fwhm < 0.5 * dtZ:
            fwFixed = True
            fwNote = ('the fit drove the IRF width below half a delay step; '
                      'it was fixed at %.3g %s' % (fwc, C['tu']))
            continue
        undet, longT = _undet_long(F, tEnd)
        if not useInf and np.count_nonzero(longT) == 1 and F.tau.size >= 2:
            Fi = auto_fit_ga(Dv, tv, F.tau[~longT], F.t0, F.fwhm, True, irf,
                             fwFixed, ca_k)
            if Fi is not None and Fi.rms <= 1.02 * F.rms:
                F = Fi
                useInf = True
                undet, longT = _undet_long(F, tEnd)
        why = ''
        if np.any(undet & ~longT):
            why = 'a time constant is not determined'
        elif F.hasInf and np.any(longT):
            why = ('a time constant beyond the time range besides the '
                   'non-decaying component')
        elif np.count_nonzero(longT) > 1:
            why = 'two time constants beyond the time range'
        elif np.any(F.tau < 0.5 * dtZ):
            why = 'a time constant below half a delay step'
        elif ca_k > 0 and np.any(F.tau < opt['caTauMin'] * F.fwhm):
            why = ('a time constant too far below the IRF width to be told '
                   'from the artifact terms')
        elif seed is not None and (F.fwhm > 3 * seed['F'].fwhm
                                   or F.fwhm < seed['F'].fwhm / 3
                                   or abs(F.t0 - seed['F'].t0) > 2 * seed['F'].fwhm):
            why = 'the IRF does not stay at time zero'
        gain = np.nan
        if n > 1:
            gain = 1 - F.rms / tabRms[n - 2]
        if n > 1 and not why and gain < opt['minGain']:
            why = 'RMS gain below %.0f %%' % (100 * opt['minGain'])
        if n == 1 or atSeed:
            why = ''
        put(n, F, gain, why)
        if not decided:
            if not why:
                nStar = n
            else:
                stopWhy = why
                decided = True
        if decided and not opt['tryAll']:
            break
        prevTau, t0c, fwc = F.tau, F.t0, F.fwhm
        n += 1
    return Result(N=tabN, Rms=tabRms, Gain=tabGain, Inf=tabInf, Tau=tabTau,
                  Std=tabStd, Why=tabWhy, F=tabF, nStar=nStar, stopWhy=stopWhy,
                  fwFixed=fwFixed, fwNote=fwNote, caK=int(ca_k))


def auto_oscillations(wl, t, R, t0, fwhm, sgw=None, time_unit='ps'):
    """Oscillations left in the residual R (n_wavelength x n_delay) of the
    chosen fit: noise-weighted power spectrum after the IRF, its peaks
    and, for every peak, the wavelengths of largest amplitude and the
    node between them. pkIdx holds 0-based indices into freq."""
    from .coherence import compute_coherence
    O = Result(status='none', msg='', win=[np.nan, np.nan], res=np.nan,
               fMax=np.nan, freq=np.zeros(0), power=np.zeros(0),
               floor=np.zeros(0), pkF=np.zeros(0), pkRatio=np.zeros(0),
               pkIdx=np.zeros(0, int), pkWl=[], pkNode=np.zeros(0), Pmap=None)
    if time_unit not in ('fs', 'ps'):
        O.status = 'skipped'
        O.msg = 'the time unit is not fs or ps'
        return O
    if t0 is None or not np.isfinite(t0):
        O.status = 'skipped'
        O.msg = 'no time zero in the data'
        return O
    wl = np.asarray(wl, float).ravel()
    t = np.asarray(t, float).ravel()
    R = np.asarray(R, float)
    hit = np.flatnonzero(t >= t0 + 2.5 * fwhm)
    if hit.size == 0 or hit[0] + 1 > t.size - 32:
        O.status = 'short'
        return O
    jA = int(hit[0]) + 1                       # 1-based
    d = np.diff(t[jA - 1:])
    d0 = np.median(d[:min(5, d.size)])
    k = 1
    while k <= d.size and d[k - 1] <= 1.5 * d0:
        k += 1
    jB = jA + k - 1
    if jB - jA + 1 < 32:
        O.status = 'short'
        return O
    O.win = [float(t[jA - 1]), float(t[jB - 1])]
    try:
        ci = compute_coherence(R, t, t[jA - 1], t[jB - 1], 'hann', 4, 'linear',
                               'cm-1', 'none', time_unit)
    except Exception as errO:
        O.status = 'failed'
        O.msg = str(errO)
        return O
    Pm = np.asarray(ci.pmap, float)
    fr = np.asarray(ci.freq, float).ravel()
    wgt = np.ones(Pm.shape[0])
    if sgw is not None and np.size(sgw) == Pm.shape[0]:
        sgw = np.asarray(sgw, float).ravel()
        with np.errstate(invalid='ignore'):
            okw = np.isfinite(sgw) & (sgw > 0)
        if np.any(okw):
            wgt[okw] = (np.median(sgw[okw]) / sgw[okw]) ** 2
            wgt[~okw] = 0
    Ps = np.sum(Pm * wgt[:, None], axis=0)
    res = 2 * fr[-1] / ci.n_win
    df = fr[1] - fr[0]
    hw = int(max(10, matlab_round(10 * res / df)))
    flr = _run_med(Ps, 2 * hw + 1)
    nF = fr.size
    isMax = np.zeros(nF, bool)
    isMax[1:nF - 1] = (Ps[1:nF - 1] >= Ps[:nF - 2]) & (Ps[1:nF - 1] > Ps[2:nF])
    inR = fr >= 2 * res
    dev = Ps[inR] - flr[inR]
    spr = _mad(dev)
    with np.errstate(invalid='ignore'):
        cand = np.flatnonzero(isMax & inR & (Ps > flr + 10 * spr) & (Ps > 2 * flr))
    cand = cand[np.argsort(-Ps[cand], kind='stable')]
    keep = []
    for c in cand:
        if not keep or np.all(np.abs(fr[c] - fr[keep]) >= 2 * res):
            keep.append(int(c))
        if len(keep) >= 6:
            break
    keep = np.array(sorted(keep), int)
    O.res, O.fMax = float(res), float(fr[-1])
    O.freq, O.power, O.floor, O.Pmap = fr, Ps, flr, Pm
    O.pkIdx = keep
    O.pkF = fr[keep]
    with np.errstate(invalid='ignore', divide='ignore'):
        O.pkRatio = Ps[keep] / flr[keep]
    O.pkWl = []
    O.pkNode = np.full(keep.size, np.nan)
    for q, kq in enumerate(keep):
        pr = Pm[:, kq]
        pr = (pr + np.concatenate([[pr[0]], pr[:-1]])
              + np.concatenate([pr[1:], [pr[-1]]])) / 3
        i1 = int(np.argmax(pr))
        m1 = pr[i1]
        bestJ, bestV, bestNode = -1, 0.0, np.nan
        for j in range(1, pr.size - 1):
            if j == i1 or pr[j] < pr[j - 1] or pr[j] < pr[j + 1]:
                continue
            if pr[j] < 0.25 * m1:
                continue
            a, b = min(i1, j), max(i1, j)
            seg = pr[a:b + 1]
            idd = int(np.argmin(seg))
            if seg[idd] < 0.3 * pr[j] and pr[j] > bestV:
                bestJ, bestV, bestNode = j, pr[j], wl[a + idd]
        if bestJ >= 0:
            O.pkWl.append(np.sort(np.array([wl[i1], wl[bestJ]])))
            O.pkNode[q] = bestNode
        else:
            O.pkWl.append(np.array([wl[i1]]))
    O.status = 'nopeak' if keep.size == 0 else 'found'
    return O


def auto_select_model(P, S, n):
    """Make the fit with n components the model of the summary S: spectra,
    residual, how alike consecutive spectra are, oscillations."""
    from .sequential import compute_eads_from_dads
    n = int(max(1, min(len(S.tabF), matlab_round(n))))
    F = S.tabF[n - 1]
    wl = np.asarray(P.wl, float).ravel()
    rows, tv, wv = S.rows, S.t, S.wl
    Dv = P.D[np.ix_(rows, S.cols)]
    tEnd = S.tEnd
    undet, longT = _undet_long(F, tEnd)
    longNote = ''
    if F.hasInf:
        longNote = ('one component does not decay within the time range; '
                    'its lifetime cannot be given')
    elif np.any(longT):
        kL = int(np.flatnonzero(longT)[0])
        if undet[kL]:
            longNote = ('the time constant %.3g %s lies beyond the time range '
                        'and is not determined' % (F.tau[kL], P.tu))
        else:
            longNote = ('the time constant %.3g %s is longer than twice the '
                        'time range: treat it with caution' % (F.tau[kL], P.tu))
    nc = F.tau.size + int(F.hasInf)
    DADS = np.full((wl.size, nc), np.nan)
    DADS[rows, :] = F.A
    EADS = np.full((wl.size, nc), np.nan)
    try:
        EADS[rows, :] = compute_eads_from_dads(F.A, F.tau, F.hasInf)[0]
    except Exception:
        pass
    caSpec = np.full((wl.size, F.caK), np.nan)
    caSpec[rows, :] = F.caSpec
    caPkW = np.full(F.caK, np.nan)
    caPkV = np.full(F.caK, np.nan)
    for k in range(F.caK):
        ii = int(np.argmax(np.abs(F.caSpec[:, k])))
        caPkW[k], caPkV[k] = wv[ii], F.caSpec[ii, k]
    Rv = Dv - F.fit
    Ur, _, Vrt = np.linalg.svd(Rv, full_matrices=False)
    simR = np.full(max(0, nc - 1), np.nan)
    simS = simR.copy()
    for k in range(nc - 1):
        e1, e2 = EADS[rows, k], EADS[rows, k + 1]
        if np.all(np.isfinite(e1)) and np.all(np.isfinite(e2)) and \
                np.linalg.norm(e1) > 0 and np.linalg.norm(e2) > 0:
            simR[k] = float(e1 @ e2) / (np.linalg.norm(e1) * np.linalg.norm(e2))
            simS[k] = float(e1 @ e2) / float(e1 @ e1)
    pkW = np.full(nc, np.nan)
    pkV = np.full(nc, np.nan)
    for k in range(nc):
        ii = int(np.argmax(np.abs(F.A[:, k])))
        pkW[k], pkV[k] = wv[ii], F.A[ii, k]
    S.model = Result(tau=F.tau, tauStd=F.tauStd, t0=F.t0, t0Std=F.t0Std,
                     fwhm=F.fwhm, fwhmStd=F.fwhmStd, hasInf=F.hasInf,
                     irf=S.irf, fwFixed=S.fwFixed, irfSeed=F.irfFixed,
                     DADS=DADS, EADS=EADS, rms=F.rms, fit=F.fit, resid=Rv,
                     acT=auto_ac1(Vrt[0, :]), acW=auto_ac1(Ur[:, 0]),
                     simR=simR, simS=simS, pkW=pkW, pkV=pkV,
                     longNote=longNote, caK=F.caK, caSpec=caSpec,
                     caPkW=caPkW, caPkV=caPkV)
    S.nUse = n
    sgw = np.asarray(P.noiseOut, float)[rows]
    if S.irf == 'numerical':
        S.osc = auto_oscillations(wv, tv, Rv, F.t0, F.fwhm, sgw, P.tu)
    else:
        S.osc = auto_oscillations(wv, tv, Rv, np.nan, np.nan, sgw, P.tu)
    return S


def auto_summarize(P, opt=None, prog=None):
    """Singular values, global fits with a growing number of components
    (with coherent-artifact terms if they are needed), the chosen model.
    prog(fraction, message) is called along the way."""
    opt = auto_options(opt)
    wl = np.asarray(P.wl, float).ravel()
    t = np.asarray(P.t, float).ravel()
    D = np.asarray(P.D, float)
    rows = np.any(np.isfinite(D), axis=1)
    cols = np.all(np.isfinite(D[rows, :]), axis=0)
    if np.count_nonzero(rows) < 5 or \
            np.count_nonzero(cols) < max(10, 0.5 * t.size):
        raise ValueError('Too many missing values for the automatic '
                         'analysis. Crop the data to the region that is '
                         'filled completely.')
    Dv = D[np.ix_(rows, cols)]
    tv = t[cols]
    wv = wl[rows]
    sgv = np.asarray(P.noiseOut, float)[rows]
    S = Result()
    S.rows, S.cols, S.t, S.wl = rows, cols, tv, wv
    S.rmsNoise = float(np.sqrt(np.mean(sgv[np.isfinite(sgv)] ** 2)))
    U, sv, Vt = np.linalg.svd(Dv, full_matrices=False)
    nsv = sv.size
    fl = float(np.median(sv[math.ceil(nsv / 2) - 1:]))
    if not (fl > 0):
        fl = float(np.finfo(float).eps)
    kMax = min(10, nsv)
    acU = np.array([auto_ac1(U[:, k]) for k in range(kMax)])
    acV = np.array([auto_ac1(Vt[k, :]) for k in range(kMax)])
    nStruct = 0
    while nStruct < kMax and acU[nStruct] > 0.5 and acV[nStruct] > 0.5:
        nStruct += 1
    S.sv, S.svFloor = sv, fl
    S.nAbove = int(np.count_nonzero(sv > 5 * fl))
    S.acU, S.acV, S.nStruct = acU, acV, nStruct
    dAll = np.diff(tv)
    dtMin = float(np.min(dAll[dAll > 0]))
    dtZ = P.chirp['dtZ']
    hasT0 = bool(P.hasT0 and np.isfinite(P.t0Guess)
                 and np.count_nonzero(tv < P.t0Guess) >= 3)
    if hasT0:
        irf = 'numerical'
        t0s = float(P.t0Guess)
        if not np.isfinite(dtZ) or dtZ <= 0:
            dtZ = dtMin
        fws = P.fwGuess
        if not np.isfinite(fws) or fws <= 0:
            fws = 2 * dtZ
        fws = max(fws, dtZ)
        tauLo = max(0.5 * fws, dtZ)
    else:
        irf = 'none'
        t0s = min(0.0, float(tv[0]))
        dtZ = dtMin
        fws = dtMin
        tauLo = 2 * dtMin
    tEnd = float(tv[-1] - t0s)
    tauHi = max(0.5 * tEnd, 8 * tauLo)
    S.irf, S.tEnd, S.dtZ = irf, tEnd, dtZ
    kc = min(10, nsv)
    Dk = sv[:kc, None] * Vt[:kc, :]
    C = dict(irf=irf, t0s=t0s, fws=fws, dtZ=dtZ, tauLo=tauLo, tauHi=tauHi,
             tEnd=tEnd, tu=P.tu)
    T = auto_fit_sequence(Dv, Dk, tv, C, opt, 0, prog, 0.15, 0.50)
    if T.nStar == 0:
        raise ValueError('The global fits of the automatic analysis failed.')
    ca = Result(mode=opt['ca'], use=0, why='', gain=np.nan, note='',
                tauFast=np.nan, fw0=np.nan, n0=T.nStar,
                rms0=T.Rms[T.nStar - 1], nSeed=0)
    nPre = 0
    if irf != 'none':
        ca.fw0 = T.F[T.nStar - 1].fwhm
        nPre = int(np.count_nonzero(tv < T.F[T.nStar - 1].t0 - 1.5 * ca.fw0))
    if irf == 'none':
        ca.note = 'not applicable: no time zero in the data'
    elif opt['ca'] == 'auto' and nPre < 5:
        ca.note = ('not tried: only %d delay(s) before the rise, so a baseline '
                   'error could not be told from an artifact' % nPre)
    elif opt['ca'] != 'off':
        if prog is not None:
            prog(0.52, 'Coherent artifact...')
        F0 = T.F[T.nStar - 1]
        fast = F0.tau < opt['caTauMin'] * F0.fwhm
        tauK = F0.tau[~fast]
        if tauK.size == 0:
            tauK = F0.tau
        F3 = None
        for _ in range(3):
            F3 = auto_fit_ga(Dv, tv, tauK, F0.t0, F0.fwhm, F0.hasInf, irf,
                             T.fwFixed, 3)
            if F3 is None or F3.tau.size < 2:
                break
            with np.errstate(invalid='ignore', divide='ignore'):
                rel1 = F3.tauStd[0] / F3.tau[0]
            if F3.tau[0] >= opt['caTauMin'] * F3.fwhm and np.isfinite(rel1) \
                    and rel1 <= 0.5:
                break
            tauK = F3.tau[1:]
        if F3 is None:
            ca.note = ('the fit with artifact terms failed (too few delays '
                       'around time zero?)')
        elif F3.t0 < tv[0] or F3.t0 > tv[-1] or F3.fwhm > 3 * F0.fwhm or \
                F3.fwhm < F0.fwhm / 3:
            ca.note = ('with artifact terms the fitted IRF does not stay at '
                       'time zero')
        else:
            ca.gain = 1 - F3.rms / F0.rms
            why = ''
            if opt['ca'] == 'on':
                why = 'on'
            elif np.any(fast):
                why = 'fast'
                ca.tauFast = float(np.min(F0.tau[fast]))
            elif ca.gain >= opt['caGain']:
                why = 'gain'
            if why:
                C['seed'] = dict(n=F3.tau.size + int(F3.hasInf), F=F3,
                                 tau0=T.Tau, inf0=T.Inf, fwFixed=T.fwFixed,
                                 fwNote=T.fwNote)
                T3 = auto_fit_sequence(Dv, Dk, tv, C, opt, 3, prog, 0.55, 0.85)
                if T3.nStar > 0:
                    T = T3
                    ca.use = 3
                    ca.why = why
                    ca.nSeed = C['seed']['n']
                else:
                    ca.note = 'the fits with artifact terms failed'
    S.ca = ca
    S.tabN, S.tabRms, S.tabGain = T.N, T.Rms, T.Gain
    S.tabInf, S.tabTau, S.tabStd = T.Inf, T.Tau, T.Std
    S.tabWhy, S.tabF = T.Why, T.F
    S.nStar, S.fwNote, S.fwFixed = T.nStar, T.fwNote, T.fwFixed
    S.stopWhy = T.stopWhy
    S.hitMax = not T.stopWhy
    S.minGain = opt['minGain']
    S.caGainMin = opt['caGain']
    if prog is not None:
        prog(0.9, 'Chosen model...')
    S = auto_select_model(P, S, S.nStar)
    if prog is not None:
        prog(1, 'Done.')
    return S


# ======================================================================
# report
# ======================================================================
_NONFIN = None


def _sp(fmt, *args):
    """sprintf with MATLAB's spelling of NaN and Inf."""
    global _NONFIN
    s = fmt % args
    for a in args:
        if isinstance(a, (float, np.floating)) and not np.isfinite(a):
            if _NONFIN is None:
                import re
                _NONFIN = (re.compile(r'\bnan\b'), re.compile(r'\binf\b'))
            s = _NONFIN[0].sub('NaN', s)
            s = _NONFIN[1].sub('Inf', s)
            break
    return s


def auto_report_lines(P, S=None, src='current data', version='1.0.2'):
    """Text of the report of the automatic analysis (list of lines):
    corrections; with a summary S also the number of components, the
    chosen fit, oscillations and the list of things to check by eye."""
    if not src:
        src = 'current data'
    tu, xu = P.tu, P.xu
    wlIn = np.asarray(P.wlIn, float).ravel()
    tIn = np.asarray(P.tIn, float).ravel()
    L = []
    L.append(_sp('AUTOMATIC ANALYSIS   (TRSpecAnalyzer v%s)', version))
    L.append(_sp('Data: %s', src))
    L.append(_sp('      %d wavelengths (%.4g ... %.4g %s) x %d delays '
                 '(%.4g ... %.4g %s)', wlIn.size, np.min(wlIn), np.max(wlIn),
                 xu, tIn.size, np.min(tIn), np.max(tIn), tu))
    L.append('')
    L.append('1. CORRECTIONS (proposal)')
    ok = np.asarray(P.rowOK, bool)
    nz = np.asarray(P.noise, float)[ok]
    L.append(_sp('  Noise        median %.3g, from %.3g to %.3g along the '
                 'spectrum', _median(nz), np.min(nz), np.max(nz)))
    if P.hasT0:
        L.append(_sp('  Signal       starts at t = %.4g %s (delay no. %d)',
                     P.tOn, tu, int(P.jOn)))
    else:
        L.append(_sp('  Signal       no time zero found: %s', P.t0Why))
    if P.bgStatus == 'applied':
        ab = np.abs(np.asarray(P.bg, float))
        ib = int(np.nanargmax(ab))
        L.append(_sp('  Background   mean of the first %d delays '
                     '(t <= %.4g %s) subtracted', P.nBg, tIn[P.nBg - 1], tu))
        L.append(_sp('               largest value %.3g at %.4g %s; '
                     'on average %.1f times its own uncertainty',
                     ab[ib], wlIn[ib], xu, P.bgSig))
    elif P.bgStatus == 'zero':
        L.append(_sp('  Background   none: the mean of the first %d '
                     'delays is zero within the noise', P.nBg))
    elif P.bgStatus == 'few':
        L.append('  Background   not set: too few delays before the signal')
    elif P.bgStatus == 'signal':
        L.append('  Background   not set: the first delays look like signal')
    elif P.bgStatus == 'notime0':
        L.append('  Background   not set: no time zero')
    else:
        L.append('  Background   switched off')
    maskReg = np.asarray(P.maskReg, float).reshape(-1, 4)
    keptReg = np.asarray(P.keptReg, float).reshape(-1, 4)
    if P.maskStatus == 'applied':
        if maskReg.shape[0] == 0:
            L.append('  Mask         none inside the spectrum (noisy end: '
                     'see Crop)')
        for q in range(maskReg.shape[0]):
            if np.isfinite(P.opt['maskSN']):
                L.append(_sp('  Mask         %.4g ... %.4g %s '
                             '(noise up to %.1f times the neighbors, signal '
                             'down to %.0f times the noise)',
                             maskReg[q, 0], maskReg[q, 1], xu, maskReg[q, 2],
                             maskReg[q, 3]))
            else:
                L.append(_sp('  Mask         %.4g ... %.4g %s '
                             '(noise up to %.1f times the neighbors)',
                             maskReg[q, 0], maskReg[q, 1], xu, maskReg[q, 2]))
    elif P.maskStatus == 'none':
        if keptReg.shape[0] == 0:
            L.append(_sp('  Mask         none: no channel is more '
                         'than %.3g times noisier than its neighbors',
                         P.opt['maskK']))
        else:
            L.append('  Mask         none: the noisy channels show a '
                     'clear signal')
    elif P.maskStatus == 'toomany':
        L.append('  Mask         not set: too many noisy channels')
    else:
        L.append('  Mask         switched off')
    if keptReg.shape[0]:
        L.append('               noisier than the neighbors, but kept:')
        nK = keptReg.shape[0]
        for q in range(min(nK, 6)):
            L.append(_sp('                 %.4g ... %.4g %s  (noise up to '
                         '%.1f times, signal at least %.0f times the noise)',
                         keptReg[q, 0], keptReg[q, 1], xu, keptReg[q, 2],
                         keptReg[q, 3]))
        if nK > 6:
            L.append(_sp('                 ... and %d more ranges', nK - 6))
    ch = P.chirp
    vr = ok & ~np.asarray(P.maskRows, bool)
    if ch['status'] in ('applied', 'found'):
        wv = wlIn[vr]
        tv0 = np.asarray(ch['t0'], float)[vr]
        L.append(_sp('  Chirp        time zero runs from %.3g %s at '
                     '%.4g %s to %.3g %s at %.4g %s (span %.3g %s)',
                     tv0[0], tu, wv[0], xu, tv0[-1], tu, wv[-1], xu,
                     ch['span'], tu))
        if ch['model'] == 'parabola':
            L.append('               (parabola: the chirp formula gave no '
                     'sensible curve)')
        L.append(_sp('               curve through the rise of %d wavelengths '
                     '(RMS %.2g %s); %.0f %% of all rises within %.2g %s of it',
                     np.asarray(ch['pts']).reshape(-1, 2).shape[0], ch['rms'],
                     tu, 100 * ch['fracIn'], ch['tol'], tu))
        if ch['status'] == 'applied':
            L.append('               corrected: time zero is now 0 at every '
                     'wavelength')
        else:
            L.append('               NOT corrected (switched off)')
    elif ch['status'] == 'none':
        L.append(_sp('  Chirp        none: time zero varies by %.2g %s along '
                     'the spectrum (negligible below %.2g %s)',
                     ch['span'], tu, 0.75 * ch['tol'], tu))
        if ch['shift'] != 0:
            L.append(_sp('               time axis shifted by %.4g %s so that '
                         'time zero is 0', -ch['shift'], tu))
    elif ch['status'] == 'failed':
        L.append(_sp('  Chirp        not determined: %s', ch['msg']))
    elif ch['status'] == 'notime0':
        L.append('  Chirp        not determined: no time zero')
    else:
        L.append('  Chirp        switched off')
    if P.cropT is None and P.cropW is None:
        if P.opt['doCrop']:
            L.append('  Crop         nothing to remove')
        else:
            L.append('  Crop         switched off')
    else:
        if P.cropT is not None:
            L.append(_sp('  Crop         t = %.4g ... %.4g %s (%d delays '
                         'without data at all wavelengths removed)',
                         P.cropT[0], P.cropT[1], tu,
                         int(np.count_nonzero(~np.asarray(P.keepT, bool)))))
        if P.cropW is not None:
            L.append(_sp('  Crop         %.4g ... %.4g %s (%d noisy channels '
                         'at the end of the spectrum removed)',
                         P.cropW[0], P.cropW[1], xu,
                         int(np.count_nonzero(~np.asarray(P.keepW, bool)))))
    L.append(_sp('  Result       %d wavelengths x %d delays',
                 np.size(P.wl), np.size(P.t)))
    if S is None:
        return L

    L.append('')
    L.append('2. NUMBER OF COMPONENTS')
    kS = min(8, S.sv.size)
    L.append(_sp('  Singular values / first:  %s',
                 '  '.join(_sp('%.2g', S.sv[k] / S.sv[0]) for k in range(kS))))
    L.append(_sp('  above 5 x noise floor: %d;  with smooth vectors '
                 '(spectrum and time): %d', S.nAbove, S.nStruct))
    ca = S.ca
    if ca.mode == 'off':
        L.append('  Coherent artifact terms: switched off')
    elif ca.note:
        L.append(_sp('  Coherent artifact terms: none (%s)', ca.note))
    elif ca.use > 0 and ca.why == 'on':
        L.append(_sp('  Coherent artifact terms: %d, switched on '
                     '(RMS without them: %.4g with %d component(s))',
                     ca.use, ca.rms0, ca.n0))
    elif ca.use > 0 and ca.why == 'fast':
        L.append(_sp('  Coherent artifact terms: %d (Gaussian of the IRF '
                     'and its two derivatives).', ca.use))
        L.append(_sp('    Without them the fits need a time constant of '
                     '%.3g %s, below a third of the IRF width', ca.tauFast, tu))
        L.append(_sp('    (%.3g %s). A signal that short cannot be told '
                     'from a coherent artifact;', ca.fw0, tu))
        L.append('    it is described by the terms and not counted as a '
                 'component.')
        L.append(_sp('    (RMS without the terms: %.4g with %d '
                     'component(s).)', ca.rms0, ca.n0))
    elif ca.use > 0:
        L.append(_sp('  Coherent artifact terms: %d (Gaussian of the IRF '
                     'and its two derivatives). They lower the RMS by %.1f %%;',
                     ca.use, 100 * ca.gain))
        L.append('    the components below were counted with them in place.')
    else:
        L.append(_sp('  Coherent artifact terms: not needed (three terms '
                     'lower the RMS by %.1f %%, less than %.0f %%)',
                     100 * ca.gain, 100 * S.caGainMin))
    if ca.use > 0 and np.isfinite(ca.fw0):
        L.append(_sp('    IRF FWHM: %.4g %s without the terms, %.4g %s '
                     'with them.', ca.fw0, tu, S.tabF[S.nStar - 1].fwhm, tu))
    L.append('  Global fits with N components (time constants and IRF free):')
    L.append(_sp('    N   RMS        gain    time constants (%s)', tu))
    for n in range(1, len(S.tabN) + 1):
        tp = [_sp('%.3g', v) for v in np.atleast_1d(S.tabTau[n - 1])]
        if S.tabInf[n - 1]:
            tp.append('non-decaying')
        if n == 1:
            gs = '   -  '
        else:
            gs = _sp('%5.1f%%', float(100 * S.tabGain[n - 1]))
        mark = ''
        if S.tabWhy[n - 1]:
            mark = '   (' + S.tabWhy[n - 1] + ')'
        if n == S.nStar and n == S.nUse:
            mark = '   <- chosen'
        elif n == S.nStar:
            mark = '   <- automatic choice'
        elif n == S.nUse:
            mark = mark + '   <- used (set by hand)'
        L.append(_sp('    %d   %-9.4g  %s  %s%s', n, S.tabRms[n - 1], gs,
                     ', '.join(tp), mark))
    if ca.use > 0 and ca.nSeed > 1:
        L.append(_sp('  (N below %d: t0 and FWHM kept at the values of the '
                     'fit with %d, for comparison.)', ca.nSeed, ca.nSeed))
    if S.hitMax:
        L.append(_sp('  Automatic choice: %d, the largest number tried. A '
                     'further component may still improve the fit.', S.nStar))
    else:
        L.append(_sp('  Automatic choice: %d. One more was rejected: %s.',
                     S.nStar, S.stopWhy))
    if S.nUse != S.nStar:
        L.append(_sp('  Used below: %d components, set by hand '
                     '("Components used").', S.nUse))
    Mo = S.model
    L.append('')
    if S.nUse == S.nStar:
        L.append('3. CHOSEN FIT')
    else:
        L.append(_sp('3. FIT WITH %d COMPONENTS (set by hand)', S.nUse))
    sign = lambda v: 'positive' if v >= 0 else 'negative'
    for k in range(Mo.tau.size):
        if np.isfinite(Mo.tauStd[k]):
            es = _sp('+/- %.2g', Mo.tauStd[k])
        else:
            es = '[not determined]'
        L.append(_sp('  tau_%d = %-9.4g %-16s %s   DADS largest at '
                     '%.4g %s (%s)', k + 1, Mo.tau[k], es, tu, Mo.pkW[k], xu,
                     sign(Mo.pkV[k])))
    if Mo.hasInf:
        L.append(_sp('  non-decaying component              DADS largest at '
                     '%.4g %s (%s)', Mo.pkW[-1], xu, sign(Mo.pkV[-1])))
    if Mo.longNote:
        L.append(_sp('  Note: %s.', Mo.longNote))
    if Mo.irf == 'none':
        L.append('  IRF: none (pure exponentials; no time zero in the data)')
    elif Mo.irfSeed:
        L.append(_sp('  IRF t0   = %.3g %s (fixed)', Mo.t0, tu))
        L.append(_sp('  IRF FWHM = %.4g %s (fixed)', Mo.fwhm, tu))
        L.append(_sp('  Note: t0 and FWHM are those of the fit with %d '
                     'components. In a fit with too few', ca.nSeed))
        L.append('        components, free artifact terms leave time zero '
                 'and stand in for the missing kinetics.')
    else:
        if np.isfinite(Mo.t0Std):
            L.append(_sp('  IRF t0   = %.3g +/- %.2g %s', Mo.t0, Mo.t0Std, tu))
        else:
            L.append(_sp('  IRF t0   = %.3g %s', Mo.t0, tu))
        if Mo.fwFixed:
            L.append(_sp('  IRF FWHM = %.4g %s (fixed)', Mo.fwhm, tu))
            L.append(_sp('  Note: %s.', S.fwNote))
        elif np.isfinite(Mo.fwhmStd):
            L.append(_sp('  IRF FWHM = %.4g +/- %.2g %s', Mo.fwhm, Mo.fwhmStd,
                         tu))
        else:
            L.append(_sp('  IRF FWHM = %.4g %s', Mo.fwhm, tu))
    if Mo.caK > 0:
        caNm = ('Gaussian', '1st derivative', '2nd derivative')
        L.append(_sp('  Coherent artifact: %d terms of the IRF width, '
                     'spectra free', Mo.caK))
        for k in range(Mo.caK):
            L.append(_sp('    %-14s largest %.3g at %.4g %s', caNm[k],
                         Mo.caPkV[k], Mo.caPkW[k], xu))
    L.append(_sp('  RMS = %.4g   (noise estimate %.3g)', Mo.rms, S.rmsNoise))
    if Mo.acT > 0.5 and Mo.acW > 0.5:
        L.append(_sp('  Residual: structured (autocorrelation of its main '
                     'part: %.2f in time, %.2f along the spectrum).',
                     Mo.acT, Mo.acW))
        L.append('            Not everything is described by the exponentials.')
    else:
        L.append(_sp('  Residual: no clear structure (autocorrelation %.2f '
                     'in time, %.2f along the spectrum).', Mo.acT, Mo.acW))
    if Mo.simR.size >= 1:
        L.append('  EADS (sequential picture 1 -> 2 -> ...), consecutive '
                 'spectra:')
        for k in range(Mo.simR.size):
            if not np.isfinite(Mo.simR[k]):
                continue
            if Mo.simR[k] > 0.95:
                hint = 'same shape'
            elif Mo.simR[k] > 0.8:
                hint = 'similar shape'
            else:
                hint = 'different shape'
            L.append(_sp('    EADS %d -> %d: shape correlation %.2f (%s), '
                         'size ratio %.2f', k + 1, k + 2, Mo.simR[k], hint,
                         Mo.simS[k]))
    O = S.osc
    L.append('')
    L.append('4. OSCILLATIONS IN THE RESIDUAL')
    if O.status in ('found', 'nopeak'):
        L.append(_sp('  Window %.4g ... %.4g %s; detectable %.0f ... %.0f cm-1 '
                     '(resolution about %.0f cm-1)', O.win[0], O.win[1], tu,
                     2 * O.res, O.fMax, O.res))
        if O.status == 'nopeak':
            L.append('  No significant peak.')
        for q in range(O.pkF.size):
            ws = O.pkWl[q]
            if np.size(ws) == 2:
                wtxt = _sp('strongest at %.4g and %.4g %s, node near %.4g %s',
                           ws[0], ws[1], xu, O.pkNode[q], xu)
            else:
                wtxt = _sp('strongest at %.4g %s', ws[0], xu)
            L.append(_sp('  %6.0f cm-1   %5.1f x floor   %s', O.pkF[q],
                         O.pkRatio[q], wtxt))
    elif O.status == 'short':
        L.append('  Not analyzed: fewer than 32 evenly spaced delays after '
                 'time zero.')
    elif O.status == 'skipped':
        L.append(_sp('  Not analyzed: %s.', O.msg))
    else:
        L.append(_sp('  Not analyzed: %s', O.msg))
    C = list(P.checks)
    if S.nUse != S.nStar and S.tabWhy[S.nUse - 1]:
        C.append(_sp('The fit with %d components was set by hand. '
                     'The automatic choice rejected it: %s.', S.nUse,
                     S.tabWhy[S.nUse - 1]))
    if Mo.caK > 0:
        C.append('Artifact terms absorb every signal confined to time '
                 'zero: a real artifact, a leftover chirp, and processes '
                 'faster than the IRF. Time constants below the IRF width '
                 'are not reliable next to them.')
        fwA = S.tabF[S.nStar - 1].fwhm
        if np.isfinite(ca.fw0) and (fwA > 1.3 * ca.fw0 or fwA < ca.fw0 / 1.3):
            C.append(_sp('The artifact terms change the IRF FWHM from '
                         '%.3g to %.3g %s. With an artifact in the data the '
                         'second value is the better one. If the dashed '
                         'spectra on the Components tab do not look like an '
                         'artifact, run again with "Coherent artifact: '
                         'Off" and compare.', ca.fw0, fwA, tu))
        if ch['status'] in ('applied', 'found'):
            C.append('With an artifact in the data the steepest rise '
                     'is not exactly time zero: where an artifact of one '
                     'sign is as large as the signal it lies up to about '
                     '0.3 IRF widths off, and the chirp curve of this run '
                     'with it. If a scan of the solvent or substrate alone '
                     'exists, take the chirp from it (Chirp Corr window).')
    if not np.isfinite(P.jOn):
        C.append('The data do not change with time beyond the noise. '
                 'The fit below then describes noise or a constant signal.')
    if S.nStruct > S.nUse:
        if Mo.caK > 0:
            extra = ('the coherent artifact (described by the '
                     'artifact terms), oscillations')
        else:
            extra = 'oscillations'
        C.append(_sp('%d singular vectors are smooth but %d '
                     'component(s) were fitted. The extra ones can be '
                     '%s, a time zero that still varies along the '
                     'spectrum, or kinetics that are not exponential.',
                     S.nStruct, S.nUse, extra))
    if S.hitMax:
        C.append('The number of components stopped at the largest '
                 'number tried, not because a further one was rejected.')
    if Mo.fwFixed:
        C.append('The IRF width could not be fitted and was fixed. '
                 'Time constants close to it depend on that value.')
    if Mo.irf == 'numerical' and Mo.caK > 0 and np.any(Mo.tau < 2 * Mo.fwhm):
        C.append(_sp('A time constant (%.3g %s) is close to the IRF '
                     'width (%.3g %s) while artifact terms are in use. The '
                     'two are hard to separate: it may describe what the '
                     'terms leave of the artifact or of a leftover chirp. '
                     'Compare with the fit with one component less '
                     '("Components used").', np.min(Mo.tau), tu, Mo.fwhm, tu))
    elif Mo.irf == 'numerical' and np.any(Mo.tau < Mo.fwhm):
        C.append(_sp('A time constant (%.3g %s) is shorter than the '
                     'IRF width (%.3g %s): it depends on the IRF shape assumed '
                     '(Gaussian), and it can also be a coherent artifact at '
                     'time zero instead of a population.',
                     np.min(Mo.tau), tu, Mo.fwhm, tu))
    L.append('')
    L.append('5. TO CHECK BY EYE')
    if not C:
        L.append('  Nothing was flagged.')
    for c in C:
        L.append(_sp('  - %s', c))
    L.append('')
    L.append('The time constants describe the data; which states and steps')
    L.append('they belong to is not decided here (see Target Analysis).')
    return L
