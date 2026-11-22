"""Working state of the program and the calculations on it that need no
window: reading a file, the correction pipeline, labels, color limits.

Port of the property block and of the corresponding methods of
TRSpecAnalyzer.m. The property names of the MATLAB program are kept as
they are (deltaA, selWL, ...), so that the two versions can be compared
line by line.
"""
import re

import numpy as np

from trspec.chirp import chirp_model

VERSION = '1.0.2'            # version of the MATLAB program this is a port of
VERSION_DATE = '2026-11-22'
PY_VERSION = '1.0.2'         # version of this Python edition


def read_matrix(path):
    """Numeric matrix of a text file (comma, semicolon, tab or blank
    separated). Cells that hold no number become NaN, short rows are
    padded with NaN -- what MATLAB's readmatrix does for such files."""
    with open(path, 'r', encoding='utf-8-sig', errors='replace') as fh:
        lines = [ln.rstrip('\r\n') for ln in fh]
    lines = [ln for ln in lines if ln.strip() != '']
    if not lines:
        return np.zeros((0, 0))
    head = lines[0]
    if ',' in head:
        split = lambda s: s.split(',')
    elif ';' in head:
        split = lambda s: s.split(';')
    elif '\t' in head:
        split = lambda s: s.split('\t')
    else:
        split = lambda s: s.split()
    rows = [split(ln) for ln in lines]
    n = max(len(r) for r in rows)
    M = np.full((len(rows), n), np.nan)
    for i, r in enumerate(rows):
        for j, tok in enumerate(r):
            tok = tok.strip()
            if tok:
                try:
                    M[i, j] = float(tok)
                except ValueError:
                    pass
    return M


def parse_data_file(path):
    """(wl, t, A) of a data file: first row = delays, first column =
    wavelengths, the rest = signal (n_wl x n_t). Rows and columns whose
    axis value is not a number are dropped; both axes come back sorted.
    Port of parseDataFile."""
    M = read_matrix(path)
    if M.shape[0] < 3 or M.shape[1] < 3:
        raise ValueError('Data matrix is too small (must be at least 3x3).')
    t = M[0, 1:].copy()
    wl = M[1:, 0].copy()
    A = M[1:, 1:].copy()
    ok_w = np.isfinite(wl)
    ok_t = np.isfinite(t)
    if ok_w.sum() < 2 or ok_t.sum() < 2:
        raise ValueError('The first column (wavelengths) or the first row '
                         '(delays) of the file holds no numbers.')
    wl = wl[ok_w]
    t = t[ok_t]
    A = A[ok_w, :][:, ok_t]
    iw = np.argsort(wl, kind='stable')
    it = np.argsort(t, kind='stable')
    return wl[iw], t[it], A[iw, :][:, it]


def detect_axis_swap(wl, t):
    """True when the two axes of a file look interchanged. Port of
    detectAxisSwap."""
    wl = np.asarray(wl, float).ravel()
    t = np.asarray(t, float).ravel()
    wl_like_delay = bool(np.any(wl <= 0) or
                         (np.min(wl) > 0 and np.max(wl) / np.min(wl) > 1000))
    t_like_wl = bool(np.all(t > 150) and np.all(t < 2500) and
                     (np.max(t) - np.min(t) > 50))
    return wl_like_delay or t_like_wl


def _num(v):
    if np.isnan(v):
        return 'NaN'
    if np.isinf(v):
        return 'Inf' if v > 0 else '-Inf'
    return repr(float(v)) if abs(v) >= 1e16 or (v != 0 and abs(v) < 1e-5) \
        else ('%.15g' % v)


def write_data_file(path, wl, t, A, delim=','):
    """Write (wl, t, A) in the layout parse_data_file reads: corner 0,
    first row delays, first column wavelengths. Port of writeDataFile."""
    wl = np.asarray(wl, float).ravel()
    t = np.asarray(t, float).ravel()
    A = np.asarray(A, float)
    if A.shape != (wl.size, t.size):
        raise ValueError('Data size mismatch: A is %dx%d, expected %dx%d'
                         % (A.shape[0], A.shape[1], wl.size, t.size))
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(delim.join(['0'] + [_num(v) for v in t]) + '\n')
        for i in range(wl.size):
            fh.write(delim.join([_num(wl[i])] + [_num(v) for v in A[i, :]])
                     + '\n')


def write_matrix(path, M, delim=',', header=None):
    """A numeric matrix as text, one row per line (MATLAB writematrix);
    header: an optional first line."""
    M = np.atleast_2d(np.asarray(M, float))
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        if header is not None:
            fh.write(header + '\n')
        for row in M:
            fh.write(delim.join(_num(v) for v in row) + '\n')


# ----------------------------------------------------------------------
# session files: one .npz archive; arrays as arrays, the rest as JSON
# ----------------------------------------------------------------------
def _pack(obj, arrays):
    if isinstance(obj, np.ndarray):
        if obj.dtype == object:
            return {'__list__': [_pack(v, arrays) for v in obj.tolist()]}
        key = 'a%d' % len(arrays)
        arrays[key] = obj
        return {'__nd__': key}
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, complex):
        return {'__cx__': [obj.real, obj.imag]}
    if isinstance(obj, dict):
        return {'__dict__': {str(k): _pack(v, arrays) for k, v in obj.items()}}
    if isinstance(obj, tuple):
        return {'__tuple__': [_pack(v, arrays) for v in obj]}
    if isinstance(obj, list):
        return [_pack(v, arrays) for v in obj]
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    raise TypeError('Cannot store a value of type %s in a session.'
                    % type(obj).__name__)


def _unpack(obj, arrays):
    if isinstance(obj, list):
        return [_unpack(v, arrays) for v in obj]
    if isinstance(obj, dict):
        if '__nd__' in obj:
            return np.array(arrays[obj['__nd__']])
        if '__dict__' in obj:
            return {k: _unpack(v, arrays) for k, v in obj['__dict__'].items()}
        if '__tuple__' in obj:
            return tuple(_unpack(v, arrays) for v in obj['__tuple__'])
        if '__list__' in obj:
            return [_unpack(v, arrays) for v in obj['__list__']]
        if '__cx__' in obj:
            return complex(obj['__cx__'][0], obj['__cx__'][1])
    return obj


def save_tree(path, tree):
    import json
    arrays = {}
    js = json.dumps(_pack(tree, arrays))
    with open(path, 'wb') as fh:
        np.savez_compressed(fh, __json__=np.array(js), **arrays)


def load_tree(path):
    import json
    try:
        with np.load(path, allow_pickle=False) as z:
            if '__json__' not in z.files:
                raise KeyError('__json__')
            js = str(z['__json__'][()])
            arrays = {k: z[k] for k in z.files if k != '__json__'}
        return _unpack(json.loads(js), arrays)
    except OSError:
        raise
    except Exception:
        raise ValueError('This file does not look like a TRSpecAnalyzer '
                         'session of the Python version.')


_GREEK = [('\\rightarrow', '→'), ('\\lambda', 'λ'),
          ('\\sigma', 'σ'), ('\\alpha', 'α'),
          ('\\omega', 'ω'), ('\\Delta', 'Δ'),
          ('\\infty', '∞'), ('\\beta', 'β'),
          ('\\cdot', '·'), ('\\tau', 'τ'), ('\\nu', 'ν'),
          ('\\mu', 'μ'), ('\\pi', 'π')]


def ui_str(s):
    """Plain-text form of a label written with TeX markup (\\lambda,
    t_0, cm^{-1}): for buttons, labels and tooltips. Port of uiStr."""
    if isinstance(s, (list, tuple)):
        return [ui_str(x) for x in s]
    if not isinstance(s, str) or not s:
        return s
    for d in range(10):
        s = re.sub(r'([A-Za-z])_%d(?![0-9])' % d,
                   lambda m, d=d: m.group(1) + chr(8320 + d), s)
    s = s.replace('^{-1}', '-1')
    s = re.sub(r'\^2(?![0-9])', '²', s)
    s = re.sub(r'([_^])\{([^{}]*)\}', r'\1\2', s)
    for a, b in _GREEK:
        s = s.replace(a, b)
    return s


class AppState:
    """Data, corrections, selections and settings (no windows)."""

    def __init__(self):
        # ---- data (current working) ----
        self.wavelength = np.zeros(0)
        self.delay = np.zeros(0)
        self.deltaA = None            # n_wl x n_t, after all corrections
        self.deltaA_raw = None        # baseline the corrections act on
        self.originalWavelength = np.zeros(0)
        self.originalDelay = np.zeros(0)
        self.originalDeltaA = None    # as loaded ("Reset all")
        self.dataSourceDesc = ''
        # ---- corrections ----
        self.bgApplied = False
        self.bgSpectrum = None
        self.bgN = 0
        self.chirpApplied = False
        self.chirpPts = np.zeros((0, 2))
        self.chirpParams = None
        self.chirpFitRMS = np.nan
        self.chirpViewTMin = None
        self.chirpViewTMax = 2.0
        self.ridgeMethod = 'max|dA/dt|'
        self.ridgeTmin = -0.5
        self.ridgeTmax = 1.5
        self.ridgeSmoothN = 5
        self.ridgeWL = np.zeros(0)
        self.ridgeT = np.zeros(0)
        self.tZeroShift = 0.0
        self.maskedRegions = []       # rows (wl1, wl2, 'nan' | 'zero')
        # ---- selection, overlays ----
        self.selWL = None
        self.selT = None
        self.specOverlays = []
        self.kinOverlays = []
        # ---- view ----
        self.delayScaleMode = 'split'     # 'linear' | 'log' | 'split'
        self.splitThreshold = 1.0
        self.mainViewTMin = None
        self.mainViewTMax = None
        self.mapColormap = 'BWR'
        self.mapZMin = None
        self.mapZMax = None
        self.specUnit = 'nm'              # 'nm' | 'cm-1'
        self.timeUnit = 'ps'              # 'fs' | 'ps' | 'ns' | 'us'
        self.prefLastDir = ''
        self.holdDataChange = False
        # ---- single-wavelength fit settings ----
        self.kfNComp = 3
        self.kfTauInit = [1.0, 10.0, 100.0]
        self.kfTauFixed = [False, False, False]
        self.kfBetaInit = [1.0, 1.0, 1.0]
        self.kfBetaFixed = [True, True, True]
        self.kfStretchOn = [False, False, False]
        self.kfHasInf = False
        self.kfT0 = 0.0
        self.kfT0Fixed = True
        self.kfFWHM = 0.15
        self.kfFWHMFixed = True
        self.kfIRFTreat = 'numerical'
        self.kfFitTMin = None
        self.kfFitTMax = None
        # ---- global analysis settings ----
        self.gaNComp = 3
        self.gaTauInit = [1.0, 10.0, 100.0]
        self.gaTauFixed = [False, False, False]
        self.gaBetaInit = [1.0, 1.0, 1.0]
        self.gaStretchOn = [False, False, False]
        self.gaIRFTreat = 'numerical'
        self.gaCA = None
        self.gaSeqModel = 'age'
        self.gaHasInf = False
        self.gaT0 = 0.0
        self.gaT0Fixed = True
        self.gaFWHM = 0.15
        self.gaFWHMFixed = True
        self.gaFitTMin = None
        self.gaFitTMax = None
        # ---- global analysis results ----
        self.clearGaResults()
        # ---- results kept for the analysis windows ----
        self.cohSaved = None
        self.ldaSaved = None
        self.mcrSaved = None
        self.svdSaved = None
        self.kfSaved = None
        self.tfSaved = None
        self.tgSaved = None
        self.tgSet = None
        self.autoSaved = None       # Auto Analysis: last run and its options

    def clearGaResults(self):
        self.gaResultTau = None
        self.gaResultBeta = None
        self.gaResultT0 = None
        self.gaResultFWHM = None
        self.gaResultHasInf = False
        self.gaResultDADS = None
        self.gaResultEADS = None
        self.gaResultFit = None
        self.gaResultRMS = np.nan
        self.gaResultTauStd = None
        self.gaResultStretchOn = None
        self.gaResultIRF = ''
        self.gaResultTRange = None
        self.gaResultCA = None

    def clearChildResults(self):
        """Drop every stored analysis result (port of
        invalidateChildResults; svdSaved is kept, as in MATLAB)."""
        self.clearGaResults()
        self.cohSaved = None
        self.ldaSaved = None
        self.mcrSaved = None
        self.kfSaved = None
        self.tfSaved = None
        self.tgSaved = None
        self.autoSaved = None

    # ------------------------------------------------------------------
    # session
    # ------------------------------------------------------------------
    SESSION_FIELDS = [
        'wavelength', 'delay', 'deltaA', 'deltaA_raw',
        'originalWavelength', 'originalDelay', 'originalDeltaA',
        'dataSourceDesc',
        'bgApplied', 'bgSpectrum', 'bgN',
        'chirpApplied', 'chirpPts', 'chirpParams', 'chirpFitRMS',
        'chirpViewTMin', 'chirpViewTMax',
        'ridgeMethod', 'ridgeTmin', 'ridgeTmax', 'ridgeSmoothN',
        'ridgeWL', 'ridgeT', 'tZeroShift', 'maskedRegions',
        'selWL', 'selT', 'specOverlays', 'kinOverlays',
        'delayScaleMode', 'splitThreshold', 'mainViewTMin', 'mainViewTMax',
        'mapColormap', 'mapZMin', 'mapZMax', 'timeUnit', 'specUnit',
        'kfNComp', 'kfTauInit', 'kfTauFixed', 'kfBetaInit', 'kfBetaFixed',
        'kfStretchOn', 'kfHasInf', 'kfT0', 'kfT0Fixed', 'kfFWHM',
        'kfFWHMFixed', 'kfIRFTreat', 'kfFitTMin', 'kfFitTMax',
        'gaNComp', 'gaTauInit', 'gaTauFixed', 'gaBetaInit', 'gaStretchOn',
        'gaIRFTreat', 'gaSeqModel', 'gaHasInf', 'gaT0', 'gaT0Fixed',
        'gaFWHM', 'gaFWHMFixed', 'gaCA', 'gaResultCA',
        'gaResultTau', 'gaResultBeta', 'gaResultT0', 'gaResultFWHM',
        'gaResultHasInf', 'gaResultDADS', 'gaResultEADS', 'gaResultFit',
        'gaResultRMS', 'gaResultTauStd', 'gaFitTMin', 'gaFitTMax',
        'gaResultStretchOn', 'gaResultIRF', 'gaResultTRange',
        'cohSaved', 'ldaSaved', 'mcrSaved', 'svdSaved',
        'kfSaved', 'tfSaved', 'tgSet', 'tgSaved', 'autoSaved']

    def sessionFields(self):
        return list(self.SESSION_FIELDS)

    def sessionSave(self, path):
        """Write the whole working state to one file."""
        import datetime
        d = {k: getattr(self, k) for k in self.SESSION_FIELDS}
        d['meta'] = {'app': 'TRSpecAnalyzer', 'version': 1,
                     'appVersion': VERSION, 'edition': 'python',
                     'pyVersion': PY_VERSION,
                     'saved': datetime.datetime.now().strftime(
                         '%Y-%m-%d %H:%M:%S')}
        save_tree(path, d)

    @staticmethod
    def sessionLoad(path):
        """Content of a session file as a dictionary. Raises ValueError
        if the file is not a session of this program."""
        d = load_tree(path)
        meta = d.get('meta') if isinstance(d, dict) else None
        if not isinstance(meta, dict) or meta.get('app') != 'TRSpecAnalyzer':
            raise ValueError('This file does not look like a '
                             'TRSpecAnalyzer session (missing meta tag).')
        return d

    def sessionRead(self, path):
        """Read a session file into the state (no window is touched)."""
        return self.sessionApply(self.sessionLoad(path))

    def sessionApply(self, d):
        meta = d['meta']
        self.clearChildResults()
        for k in self.SESSION_FIELDS:
            if k in d:
                setattr(self, k, d[k])
        self.maskedRegions = [(float(a), float(b), str(m))
                              for a, b, m in (self.maskedRegions or [])]
        self.specOverlays = [float(v) for v in (self.specOverlays or [])]
        self.kinOverlays = [float(v) for v in (self.kinOverlays or [])]
        if self.gaSeqModel not in ('age', 't0'):
            self.gaSeqModel = 'age'
        return meta

    def dataStem(self):
        """File-name stem for exports, from the description of the
        loaded data. Port of dataStem."""
        d = self.dataSourceDesc
        if not isinstance(d, str) or not d:
            return ''
        sfx = ''
        m = re.match(r'^Avg of (\d+) files \((.*)\.\.\.\)$', d)
        if m:
            d = m.group(2)
            sfx = '_avg' + m.group(1)
        elif d.startswith('Loaded: '):
            d = d[8:]
        d = d.replace('\\', '/').split('/')[-1]
        nm = d.rsplit('.', 1)[0] if '.' in d[1:] else d
        return re.sub(r'[\\/:*?"<>|\s]+', '_', nm + sfx)

    # ------------------------------------------------------------------
    # labels
    # ------------------------------------------------------------------
    def xAxLabel(self):
        return ('Wavenumber (cm$^{-1}$)' if self.specUnit == 'cm-1'
                else 'Wavelength (nm)')

    def xUnitAx(self):
        return 'cm$^{-1}$' if self.specUnit == 'cm-1' else 'nm'

    def xUnitTxt(self):
        return 'cm-1' if self.specUnit == 'cm-1' else 'nm'

    def xSym(self):
        return 'ν' if self.specUnit == 'cm-1' else 'λ'

    xSymTxt = xSym

    def tUnitTxt(self):
        return {'us': 'μs', 'ns': 'ns', 'fs': 'fs'}.get(self.timeUnit, 'ps')

    tUnitAx = tUnitTxt

    def tUnitHdr(self):
        return {'us': 'us', 'ns': 'ns', 'fs': 'fs'}.get(self.timeUnit, 'ps')

    def xHdr(self):
        return 'wavenumber_cm-1' if self.specUnit == 'cm-1' else 'wavelength_nm'

    def gaDataNoCA(self):
        """Working data without the coherent-artifact part of the stored
        global fit (port of gaDataNoCA)."""
        from trspec.kinetics import ca_basis
        D = self.deltaA
        ca = self.gaResultCA
        if isinstance(ca, dict) and ca.get('order', 0) > 0 and \
                self.gaResultT0 is not None and \
                np.shape(ca['spec'])[0] == D.shape[0]:
            B = ca_basis(self.delay, self.gaResultT0, ca['fwhm'], ca['order'])
            sp = np.array(ca['spec'], float)
            sp[~np.isfinite(sp)] = 0
            D = D - sp @ B.T
        return D

    def scaleModeLabel(self):
        return {'linear': 'Linear', 'log': 'Log'}.get(
            str(self.delayScaleMode).lower(), 'Split')

    def defaultChirpTMax(self):
        """Default upper end of the chirp window for the current time
        unit. Port of defaultChirpTMax."""
        return {'fs': 500.0, 'ns': 0.002, 'us': 2e-6}.get(self.timeUnit, 2.0)

    # ------------------------------------------------------------------
    # data
    # ------------------------------------------------------------------
    def hasData(self):
        return self.deltaA is not None and self.deltaA.size > 0

    def isCropped(self):
        return (self.originalWavelength.size > 0 and
                (self.wavelength.size != self.originalWavelength.size or
                 self.delay.size != self.originalDelay.size))

    def setLoadedState(self, wl, t, A, sourceDesc):
        """State after a new data set was loaded (the part of
        setLoadedData that touches no control)."""
        wl = np.asarray(wl, float).ravel()
        t = np.asarray(t, float).ravel()
        A = np.asarray(A, float)
        self.originalWavelength = wl.copy()
        self.originalDelay = t.copy()
        self.originalDeltaA = A.copy()
        self.dataSourceDesc = sourceDesc
        self.wavelength = wl.copy()
        self.delay = t.copy()
        self.deltaA_raw = A.copy()
        self.deltaA = A.copy()
        self.bgApplied = False
        self.bgSpectrum = None
        self.bgN = 0
        self.chirpApplied = False
        self.chirpPts = np.zeros((0, 2))
        self.chirpParams = None
        self.chirpFitRMS = np.nan
        self.chirpViewTMin = float(t[0])
        self.chirpViewTMax = float(min(self.defaultChirpTMax(), t[-1]))
        self.ridgeWL = np.zeros(0)
        self.ridgeT = np.zeros(0)
        self.tZeroShift = 0.0
        self.maskedRegions = []
        self.specOverlays = []
        self.kinOverlays = []
        # MATLAB: wl(round(end/2)), 1-based
        self.selWL = float(wl[max(int(np.floor(wl.size / 2 + 0.5)) - 1, 0)])
        idx0 = np.flatnonzero(t >= 0)
        i0 = idx0[0] if idx0.size else max(int(np.floor(t.size / 2 + 0.5)) - 1, 0)
        self.selT = float(t[i0])
        self.mainViewTMin = float(t[0])
        self.mainViewTMax = float(t[-1])
        self.mapZMin = None
        self.mapZMax = None

    def recompute(self):
        """Working data from the baseline and the corrections in force:
        background, chirp, masks. Returns True when the data changed.
        Port of recompute (the caller reacts to the change)."""
        if self.deltaA_raw is None:
            return False
        data = self.deltaA_raw.copy()
        if self.bgApplied and self.bgSpectrum is not None:
            data = data - np.asarray(self.bgSpectrum, float).reshape(-1, 1)
        if self.chirpApplied and self.chirpParams is not None:
            t = self.delay
            wl = self.wavelength
            t0 = chirp_model(self.chirpParams, wl)
            new = np.full(data.shape, np.nan)
            for i in range(wl.size):
                new[i, :] = np.interp(t + t0[i], t, data[i, :],
                                      left=np.nan, right=np.nan)
            data = new
        for wl1, wl2, mode in self.maskedRegions:
            row = (self.wavelength >= wl1) & (self.wavelength <= wl2)
            if not row.any():
                continue
            data[row, :] = np.nan if mode == 'nan' else 0.0
        changed = (self.deltaA is not None and
                   not np.array_equal(data, self.deltaA, equal_nan=True))
        self.deltaA = data
        return bool(changed)

    def getChirpBaseData(self):
        base = self.deltaA_raw
        if self.bgApplied and self.bgSpectrum is not None:
            base = base - np.asarray(self.bgSpectrum, float).reshape(-1, 1)
        return base

    # ------------------------------------------------------------------
    # view helpers
    # ------------------------------------------------------------------
    def computeZLim(self):
        """Color limits of the 2D map: the values typed in, or
        symmetric about zero at the 99.5 % quantile of |data|. Port of
        computeZLim."""
        if (self.mapZMin is not None and self.mapZMax is not None and
                self.mapZMin < self.mapZMax):
            return float(self.mapZMin), float(self.mapZMax)
        if not self.hasData():
            return -1.0, 1.0
        vals = np.abs(self.deltaA).ravel()
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            cl = 1.0
        else:
            vs = np.sort(vals)
            cl = vs[max(1, int(np.ceil(0.995 * vs.size))) - 1]
            if not np.isfinite(cl) or cl == 0:
                cl = np.max(vals)
        if not np.isfinite(cl) or cl == 0:
            cl = 1.0
        return -float(cl), float(cl)

    def clipViewToMask(self, tSub, scaleStr):
        """Axis limits for the delays tSub: the display range typed in,
        cut to the data shown (and to t > 0 on a log axis). Port of
        clipViewToMask."""
        tSub = np.asarray(tSub, float)
        tmin, tmax = self.mainViewTMin, self.mainViewTMax
        useView = tmin is not None and tmax is not None and tmin < tmax
        dataLo, dataHi = float(np.min(tSub)), float(np.max(tSub))
        if useView:
            lo, hi = max(dataLo, tmin), min(dataHi, tmax)
        else:
            lo, hi = dataLo, dataHi
        pos = tSub[tSub > 0]
        if scaleStr == 'log':
            minPos = float(np.min(pos)) if pos.size else np.finfo(float).eps
            lo = max(lo, minPos)
        if hi <= lo:
            lo, hi = dataLo, dataHi
            if scaleStr == 'log' and pos.size:
                lo = float(np.min(pos))
        if hi <= lo:
            hi = lo + np.spacing(lo) * 100
            if hi <= lo:
                hi = lo + 0.1
        return lo, hi

    def mapTitle(self, scaleStr=None):
        parts = []
        if self.bgApplied:
            parts.append('BG(N=%d)' % self.bgN)
        if self.chirpApplied:
            parts.append('chirp(Sellmeier)')
        base = ', '.join(parts) if parts else 'raw'
        if self.isCropped():
            base += ', cropped'
        if scaleStr:
            return '2D ΔA or Intensity  (%s, %s scale)' % (base, scaleStr)
        return '2D ΔA or Intensity  (%s)' % base

    def statusText(self):
        if self.wavelength.size == 0:
            return 'No data loaded'
        base = self.dataSourceDesc or 'Data'
        s = '%s   (%d λ x %d t)' % (base, self.wavelength.size,
                                        self.delay.size)
        if self.isCropped():
            s += '   [cropped]'
        return s
