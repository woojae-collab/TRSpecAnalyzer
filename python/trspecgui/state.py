"""Working state of the program and the calculations on it that need no
window: reading a file, the correction pipeline, labels, color limits.

Port of the property block and of the corresponding methods of
TRSpecAnalyzer.m. The property names of the MATLAB program are kept as
they are (deltaA, selWL, ...), so that the two versions can be compared
line by line.
"""
import copy
import re

import numpy as np

from trspec.chirp import chirp_model, chirp_shift, interp_runs

VERSION = '1.1.0'            # version of the MATLAB program this is a port of
VERSION_DATE = '2026-11-26'
PY_VERSION = '1.1.0'         # version of this Python edition


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


_CALC_CACHE = []


def cachedCalc(kind, fcn, a, b, c):
    """fcn(a, b, c), or the result kept from an earlier call of the same
    kind with the same three arguments. The corrections of the main window
    are rebuilt at every redraw (recompute); an interpolation of the whole
    data set is not repeated for that. The four last results are kept;
    they are not saved. Port of cachedCalc."""
    args = (a, b, c)

    def same(x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        return x.shape == y.shape and np.array_equal(x, y, equal_nan=True)

    for i, (k, old, out) in enumerate(_CALC_CACHE):
        if k == kind and same(old[0], a) and same(old[2], c) and same(old[1], b):
            if i:
                _CALC_CACHE.insert(0, _CALC_CACHE.pop(i))
            return out.copy()
    out = fcn(a, b, c)
    _CALC_CACHE.insert(0, (kind, tuple(np.array(x, dtype=float) for x in args),
                           np.array(out, dtype=float)))
    del _CALC_CACHE[4:]
    return out


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
        # ---- signal kind (transient absorption) ----
        self.sigKind = 'dA'           # what the working data are: 'dA' | 'dTT' | 'int'
        self.sigOrigKind = 'dA'       # what the data were when loaded
        self.sigLoadKind = 'dA'       # what loaded files contain (chosen at loading)
        self.sigAUnit = 'OD'          # unit of dA in the conversion: 'OD' | 'mOD'
        # ---- reference data set (solvent / substrate) ----
        self.refOrig = None           # as loaded, on the axes of the loaded data
        self.refDesc = ''             # its file name
        self.subApplied = False       # subtracted from the data as loaded
        self.subScale = 1.0           # factor of the reference
        self.subShift = 0.0           # its shift in time (time unit of the data)
        self.subBaked = False         # a subtraction has become part of deltaA_raw
        self.baseRebased = False      # deltaA_raw is no longer the loaded data
        self.refStage = None          # the reference at the stage of deltaA_raw
        # ---- anisotropy (two polarizations) ----
        self.aniOn = False            # perpendicular data are loaded
        self.perpOrig = None          # as loaded, on the axes of the loaded data
        self.perpDesc = ''            # their file name
        self.perpStage = None         # at the stage of deltaA_raw
        self.aniView = 'r'            # what deltaA is: 'par' | 'perp' | 'magic' | 'r'
        self.aniG = 1.0               # factor of the perpendicular data (typed)
        self.aniTail = False          # take G from tail matching instead
        self.aniTailT = [np.nan, np.nan]   # delay range of the tail matching
        self.aniThr = 3.0             # r left empty where |par + 2 G perp| < aniThr x noise
        self.aniShift = 0.0           # time shift of the perpendicular data: perp(t - aniShift)
        self.aniPerp0 = None          # the corrected perpendicular data before that shift
        self.aniPar = None            # corrected parallel data (set by recompute)
        self.aniPerp = None           # corrected perpendicular data
        self.aniGUsed = 1.0           # G in force
        self.aniStore = {}            # setup and results of the analysis
        #                               windows for the views that are not
        #                               shown (view -> fields, see setAniView)
        self.aniKeepStore = False     # a change that leaves the data of other
        #                               views alone is under way (applyAnisotropy)
        # ---- solvent coherence (Coherence window) ----
        self.solvOn = False           # take the solvent coherence out of the residual
        self.solvA = 1.0              # factor of the solvent coherence
        self.solvShift = 0.0          # its shift in time (time unit of the data)
        self.solvBand = [np.nan, np.nan]   # frequency band of the fit (unit: solvBandUnit)
        self.solvBandUnit = 'cm-1'    # 'cm-1' | 'thz'
        self.solvWL = [np.nan, np.nan]     # spectral range of the fit (NaN: all)
        self.solvMaxShift = np.nan    # largest shift the fit tries (NaN: 2 delay steps)
        # (the solvent views of the main window: the data without the
        #  solvent. They exist while the subtraction is in force in the
        #  Coherence window, and take out what that analysis subtracts.)
        self.solvView = 'data'        # 'data' | 'osc' (oscillating part, in the
        #                               time window) | 'all' (whole solvent
        #                               signal, all delays)
        self.solvData = False         # a solvent view is shown (solvView is not 'data')
        self.solvDataA = 1.0          # factor of the views
        self.solvDataShift = 0.0      # their shift in time (time unit of the data)
        self.solvDataWin = [np.nan, np.nan]   # time window of 'osc'
        self.solvDataPad = 0.0        # largest shift the coherence is continued for
        self.solvStore = {}           # setup and results of the analysis
        #                               windows kept for the solvent views
        #                               not shown
        # ---- steps that made corrected data the starting point ----
        self.corrLog = []             # one dict per step, in order
        self.corrLogOK = True         # False: steps not known (older session)
        self.pendingAuto = None       # corrections of Auto Analysis being applied
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
        self.kfBandOn = False         # Fit trace on a band of wavelengths
        self.kfBandTo = np.nan        # its other end (NaN: proposed)
        self.kfBandMode = 'integral'  # 'integral' | 'mean'
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
        self.bpSaved = None         # Band position: setup and last result

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
        invalidateChildResults; svdSaved is kept, as in MATLAB). Those
        kept for the anisotropy views that are not shown go as well,
        unless the change is one that leaves the data of other views
        alone (see applyAnisotropy)."""
        self.clearLiveResults()
        if not self.aniKeepStore:
            self.aniStore = {}
        self.solvStore = {}

    def clearLiveResults(self):
        """Drop the analysis results of the data shown."""
        self.clearGaResults()
        self.cohSaved = None
        self.ldaSaved = None
        self.mcrSaved = None
        self.kfSaved = None
        self.tfSaved = None
        self.tgSaved = None
        self.autoSaved = None
        self.bpSaved = None

    # ------------------------------------------------------------------
    # session
    # ------------------------------------------------------------------
    SESSION_FIELDS = [
        'wavelength', 'delay', 'deltaA', 'deltaA_raw',
        'originalWavelength', 'originalDelay', 'originalDeltaA',
        'dataSourceDesc',
        'sigKind', 'sigOrigKind', 'sigAUnit',
        'refOrig', 'refDesc', 'subApplied', 'subScale', 'subShift',
        'subBaked', 'baseRebased', 'corrLog', 'corrLogOK',
        'aniOn', 'perpOrig', 'perpDesc', 'aniView', 'aniG', 'aniTail',
        'aniTailT', 'aniThr', 'aniShift',
        'solvOn', 'solvA', 'solvShift', 'solvBand', 'solvBandUnit', 'solvWL',
        'solvMaxShift',
        'solvData', 'solvDataA', 'solvDataShift', 'solvDataWin', 'solvDataPad',
        'solvView', 'solvStore',
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
        'kfBandOn', 'kfBandTo', 'kfBandMode',
        'gaNComp', 'gaTauInit', 'gaTauFixed', 'gaBetaInit', 'gaStretchOn',
        'gaIRFTreat', 'gaSeqModel', 'gaHasInf', 'gaT0', 'gaT0Fixed',
        'gaFWHM', 'gaFWHMFixed', 'gaCA', 'gaResultCA',
        'gaResultTau', 'gaResultBeta', 'gaResultT0', 'gaResultFWHM',
        'gaResultHasInf', 'gaResultDADS', 'gaResultEADS', 'gaResultFit',
        'gaResultRMS', 'gaResultTauStd', 'gaFitTMin', 'gaFitTMax',
        'gaResultStretchOn', 'gaResultIRF', 'gaResultTRange',
        'cohSaved', 'ldaSaved', 'mcrSaved', 'svdSaved',
        'kfSaved', 'tfSaved', 'tgSet', 'tgSaved', 'autoSaved', 'bpSaved',
        # analysis state of the anisotropy views not shown
        'aniStore']

    # The fields that hold the setup and the results of the analysis
    # windows (all are session fields). With two polarizations this set
    # exists once per view: for the view shown in the fields themselves,
    # for the others in aniStore.
    ANALYSIS_FIELDS = [
        'kfNComp', 'kfTauInit', 'kfTauFixed', 'kfBetaInit', 'kfBetaFixed',
        'kfStretchOn', 'kfHasInf', 'kfT0', 'kfT0Fixed', 'kfFWHM',
        'kfFWHMFixed', 'kfIRFTreat', 'kfFitTMin', 'kfFitTMax',
        'kfBandOn', 'kfBandTo', 'kfBandMode',
        'gaNComp', 'gaTauInit', 'gaTauFixed', 'gaBetaInit', 'gaStretchOn',
        'gaIRFTreat', 'gaSeqModel', 'gaHasInf', 'gaT0', 'gaT0Fixed',
        'gaFWHM', 'gaFWHMFixed', 'gaCA', 'gaResultCA',
        'gaResultTau', 'gaResultBeta', 'gaResultT0', 'gaResultFWHM',
        'gaResultHasInf', 'gaResultDADS', 'gaResultEADS', 'gaResultFit',
        'gaResultRMS', 'gaResultTauStd', 'gaFitTMin', 'gaFitTMax',
        'gaResultStretchOn', 'gaResultIRF', 'gaResultTRange',
        'cohSaved', 'ldaSaved', 'mcrSaved', 'svdSaved',
        'kfSaved', 'tfSaved', 'tgSet', 'tgSaved', 'autoSaved', 'bpSaved']

    def analysisFields(self):
        return list(self.ANALYSIS_FIELDS)

    def takeAnalysisState(self):
        """Setup and results of the analysis windows as one dictionary
        (a copy)."""
        return {k: copy.deepcopy(getattr(self, k))
                for k in self.ANALYSIS_FIELDS}

    def putAnalysisState(self, S):
        """Make a dictionary of takeAnalysisState the setup and the
        results in force."""
        for k in self.ANALYSIS_FIELDS:
            if k in S:
                setattr(self, k, S[k])
        self._normAnalysisFields()

    def dropViewResults(self, views):
        """Forget what is kept for the named views (the data of these
        views have changed)."""
        for v in views:
            self.aniStore.pop(v, None)

    def closeAnalysisWindows(self):
        """(The state has no windows; the application closes them.)"""

    def _normAnalysisFields(self):
        """Types of the analysis fields that a session file may hold in
        another form."""
        self.kfBandOn = bool(self.kfBandOn)
        try:
            self.kfBandTo = float(self.kfBandTo)
        except (TypeError, ValueError):
            self.kfBandTo = np.nan
        if self.kfBandMode not in ('integral', 'mean'):
            self.kfBandMode = 'integral'
        if self.gaSeqModel not in ('age', 't0'):
            self.gaSeqModel = 'age'

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
        self.clearReference()
        self.clearAnisotropy()
        self.clearCorrLog()
        for k in self.SESSION_FIELDS:
            if k in d:
                setattr(self, k, d[k])
        if self.refOrig is not None:
            self.refOrig = np.asarray(self.refOrig, float)
            if self.refOrig.size == 0:
                self.refOrig = None
        self.refDesc = str(self.refDesc or '')
        self.subApplied = bool(self.subApplied)
        self.subBaked = bool(self.subBaked)
        self.subScale = float(self.subScale)
        self.subShift = float(self.subShift)
        if 'baseRebased' not in d:
            self.baseRebased = not (
                self.deltaA_raw is not None and self.originalDeltaA is not None
                and np.array_equal(self.deltaA_raw, self.originalDeltaA,
                                   equal_nan=True))
        self.baseRebased = bool(self.baseRebased)
        if 'corrLog' not in d:
            # older session: how the starting data came about is not
            # recorded
            self.corrLog = []
            self.corrLogOK = not self.baseRebased
        if not isinstance(self.corrLog, list):
            self.corrLog = []
        self.corrLogOK = bool(self.corrLogOK)
        self.updateRefStage()
        if self.aniView not in self.ANI_VIEWS:
            self.aniView = 'r'
        if self.perpOrig is not None:
            self.perpOrig = np.asarray(self.perpOrig, float)
            if self.perpOrig.size == 0:
                self.perpOrig = None
        self.perpDesc = str(self.perpDesc or '')
        self.aniOn = bool(self.aniOn) and self.perpOrig is not None
        self.aniG = float(self.aniG)
        self.aniTail = bool(self.aniTail)
        self.aniTailT = [float(v) for v in np.ravel(self.aniTailT)][:2]
        self.aniThr = float(self.aniThr)
        try:
            self.aniShift = float(self.aniShift)
        except (TypeError, ValueError):
            self.aniShift = 0.0
        if not np.isfinite(self.aniShift):
            self.aniShift = 0.0
        self.updatePerpStage()
        self.solvA = float(self.solvA)
        self.solvShift = float(self.solvShift)
        self.solvMaxShift = float(self.solvMaxShift)
        for k in ('solvBand', 'solvWL'):
            v = [float(x) for x in np.ravel(getattr(self, k))]
            setattr(self, k, v if len(v) == 2 else [np.nan, np.nan])
        if self.solvBandUnit not in ('cm-1', 'thz'):
            self.solvBandUnit = 'cm-1'
        try:
            win = [float(x) for x in np.ravel(self.solvDataWin)]
            self.solvDataA = float(self.solvDataA)
            self.solvDataShift = float(self.solvDataShift)
            self.solvDataPad = float(self.solvDataPad)
        except (TypeError, ValueError):
            win = []
        okWin = len(win) == 2 and bool(np.all(np.isfinite(win)))
        self.solvDataWin = win if okWin else [np.nan, np.nan]
        # The solvent views. (The first candidates of version 1.1.0 stored
        # "the coherence is out of the data" without a view, with the
        # factor of the analysis set to 0.)
        if ('solvView' not in d
                and bool(np.size(self.solvData) == 1 and self.solvData)):
            self.solvView = 'osc'
            self.solvOn = True
            self.solvA = self.solvDataA
            self.solvShift = self.solvDataShift
        if self.solvView not in self.SOLV_VIEWS:
            self.solvView = 'data'
        self.solvOn = (bool(np.size(self.solvOn) == 1 and self.solvOn)
                       and self.refOrig is not None and not self.aniOn)
        if not (self.solvOn and okWin):
            self.solvView = 'data'
        self.solvData = self.solvView != 'data'
        self._normAnalysisFields()
        # what is kept for the solvent views that are not shown
        st = self.solvStore if isinstance(self.solvStore, dict) else {}
        self.solvStore = {v: dict(st[v]) for v in self.SOLV_VIEWS
                          if self.solvOn and okWin and v != self.solvView
                          and isinstance(st.get(v), dict)}
        # what is kept for the views that are not shown
        st = self.aniStore if isinstance(self.aniStore, dict) else {}
        self.aniStore = {v: dict(st[v]) for v in self.ANI_VIEWS
                         if self.aniOn and v != self.aniView
                         and isinstance(st.get(v), dict)}
        self.maskedRegions = [(float(a), float(b), str(m))
                              for a, b, m in (self.maskedRegions or [])]
        self.specOverlays = [float(v) for v in (self.specOverlays or [])]
        self.kinOverlays = [float(v) for v in (self.kinOverlays or [])]
        # sessions written before version 1.1.0 have no signal kind:
        # their data are what the file held, taken as dA
        if 'sigKind' not in d or self.sigKind not in self.SIG_KINDS:
            self.sigKind = 'dA'
        if 'sigOrigKind' not in d or self.sigOrigKind not in self.SIG_KINDS:
            self.sigOrigKind = 'dA'
        if self.sigAUnit not in ('OD', 'mOD'):
            self.sigAUnit = 'OD'
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
        self.sigOrigKind = self.sigLoadKind
        self.sigKind = self.sigLoadKind
        self.clearReference()
        self.clearAnisotropy()
        self.clearCorrLog()
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
        data = self.pipelineBase().copy()
        if self.bgApplied and self.bgSpectrum is not None:
            data = data - np.asarray(self.bgSpectrum, float).reshape(-1, 1)
        if self.chirpApplied and self.chirpParams is not None:
            t = self.delay
            wl = self.wavelength
            t0 = chirp_model(self.chirpParams, wl)
            data = self.chirpShift(t, data, t0)
        for wl1, wl2, mode in self.maskedRegions:
            row = (self.wavelength >= wl1) & (self.wavelength <= wl2)
            if not row.any():
                continue
            data[row, :] = np.nan if mode == 'nan' else 0.0
        # a solvent view of the main window (setSolvView): what the
        # Coherence window subtracts from its residual is subtracted from
        # the data themselves, in the same time window, or the whole
        # solvent signal is, at all delays
        if self.solvData:
            data = data - self.solventInData()
        # anisotropy: the perpendicular data get the same corrections; the
        # working data are the chosen view of the two
        if (self.aniOn and self.perpStage is not None
                and self.perpStage.shape == data.shape):
            perp = self.pipeOn(self.perpStage, self.pipeState(), self.delay)
            # (a time zero that differs between the two measurements)
            self.aniPerp0 = perp
            perp = self.shiftRowsKept(perp, self.delay, self.aniShift)
            self.aniPar = data
            self.aniPerp = perp
            data, self.aniGUsed, _ = self.aniCombine(
                data, perp, self.delay, self.aniView, self.aniG,
                self.aniTail, self.aniTailT, self.aniThr)
        else:
            self.aniPar = None
            self.aniPerp = None
            self.aniPerp0 = None
        changed = (self.deltaA is not None and
                   not np.array_equal(data, self.deltaA, equal_nan=True))
        self.deltaA = data
        return bool(changed)

    # ------------------------------------------------------------------
    # reference data set: data - scale * reference(t - shift)
    # ------------------------------------------------------------------
    def clearReference(self):
        """No reference data set (new data, session about to be read)."""
        self.refOrig = None
        self.refDesc = ''
        self.refStage = None
        self.subApplied = False
        self.subScale = 1.0
        self.subShift = 0.0
        self.subBaked = False
        self.baseRebased = False
        self.solvOn = False
        self.solvA = 1.0
        self.solvShift = 0.0
        self.solvBand = [np.nan, np.nan]
        self.solvBandUnit = 'cm-1'
        self.solvWL = [np.nan, np.nan]
        self.solvMaxShift = np.nan
        self.solvData = False
        self.solvDataA = 1.0
        self.solvDataShift = 0.0
        self.solvDataWin = [np.nan, np.nan]
        self.solvDataPad = 0.0
        self.solvView = 'data'
        self.solvStore = {}

    def markRebased(self):
        """deltaA_raw was replaced by corrected data (crop, signal
        conversion, corrections of Auto Analysis): a subtracted reference
        is part of it from now on."""
        if self.subApplied:
            self.subBaked = True
        self.subApplied = False
        self.baseRebased = True
        # (a solvent view that is shown is part of the data from now on,
        # like the corrections: nothing is left to subtract, and there are
        # no views any more)
        if self.solvData:
            self.solvOn = False
            self.solvA = 0.0
            self.solvDataWin = [np.nan, np.nan]
        self.solvData = False
        self.solvView = 'data'
        self.solvStore = {}
        self.updateRefStage()
        self.updatePerpStage()

    # ------------------------------------------------------------------
    # anisotropy: parallel (the loaded data) and perpendicular data
    # ------------------------------------------------------------------
    ANI_VIEWS = ('par', 'perp', 'magic', 'r')
    ANI_VIEW_NAMES = {'par': 'Parallel', 'perp': 'Perpendicular',
                      'magic': 'Magic angle', 'r': 'Anisotropy r'}

    def clearAnisotropy(self):
        """No perpendicular data (new data, session about to be read)."""
        self.aniOn = False
        self.perpOrig = None
        self.perpDesc = ''
        self.perpStage = None
        self.aniPar = None
        self.aniPerp = None
        self.aniPerp0 = None
        self.aniGUsed = 1.0
        self.aniView = 'r'
        self.aniG = 1.0
        self.aniTail = False
        self.aniTailT = [np.nan, np.nan]
        self.aniThr = 3.0
        self.aniShift = 0.0
        self.aniStore = {}

    def updatePerpStage(self):
        """The perpendicular data at the stage of deltaA_raw (None when
        they cannot be brought there)."""
        self.perpStage = None
        if self.perpOrig is None:
            return
        X, ok = self.replayCorrections(self.perpOrig)
        if ok and self.deltaA_raw is not None and X.shape == self.deltaA_raw.shape:
            self.perpStage = X

    def setPerpendicular(self, wl, t, A, desc):
        """Take a parsed file as the perpendicular data. Its axes must be
        those of the loaded (parallel) data, as parsed or with the two
        axes interchanged. Returns (ok, msg)."""
        wl0 = self.originalWavelength
        t0 = self.originalDelay
        wl = np.asarray(wl, float).ravel()
        t = np.asarray(t, float).ravel()
        A = np.asarray(A, float)
        if not (self.sameAxis(wl, wl0) and self.sameAxis(t, t0)):
            wlT, tT, AT = t.copy(), wl.copy(), A.T.copy()
            iw = np.argsort(wlT, kind='stable')
            it = np.argsort(tT, kind='stable')
            wlT, tT, AT = wlT[iw], tT[it], AT[iw, :][:, it]
            if self.sameAxis(wlT, wl0) and self.sameAxis(tT, t0):
                A = AT
            else:
                return False, (
                    'The axes of this file differ from those of the loaded '
                    'data (%d x %d points here, %d x %d in the data). The '
                    'two polarizations must be measured at the same '
                    'wavelengths and delays.'
                    % (wl.size, t.size, wl0.size, t0.size))
        if not self.corrLogOK:
            return False, (
                'This session was saved by an earlier version, after a crop: '
                'the steps made since loading are not recorded. Press "Reset '
                'all" in the main window, then load the perpendicular data.')
        if self.solvData:
            return False, (
                'A view without the solvent is shown ("View" of the main '
                'window). Choose the view "Data" first, then load the '
                'perpendicular data.')
        self.perpOrig = A.copy()
        self.perpDesc = desc
        self.updatePerpStage()
        if self.perpStage is None:
            self.perpOrig = None
            self.perpDesc = ''
            return False, ('The perpendicular data could not be brought to '
                           'the state of the data shown.')
        self.solvOn = False
        self.clearSolvDef()
        if not self.aniOn:
            # the first perpendicular data: what the analysis windows hold
            # belongs to the parallel data, and is kept as that of the
            # view "Parallel"
            self.closeAnalysisWindows()
            self.aniStore = {}
            if self.aniView != 'par':
                self.aniStore['par'] = self.takeAnalysisState()
        else:
            # other perpendicular data: every view but the parallel one
            # changes
            self.dropViewResults(('perp', 'magic', 'r'))
        return True, ''

    def aniTailProposal(self):
        """Delay range proposed for the tail matching: the last decade
        of the delays, from the first delay that is at least one tenth of
        the last one. Only a proposal: the range has to lie where the
        anisotropy has decayed. Without a positive last delay, or with
        fewer than three delays there: the last fifth of the delay
        points."""
        t = np.asarray(self.delay, float).ravel()
        tl = [float(t[max(1, int(round(0.8 * t.size))) - 1]), float(t[-1])]
        if t[-1] > 0:
            j = np.flatnonzero(t >= t[-1] / 10)
            if j.size >= 3:
                tl = [float(t[j[0]]), float(t[-1])]
        return tl

    def aniPerpAt(self, d):
        """The corrected perpendicular data shifted by d along the
        delays: perp(t - d). (d = the shift in use: the data in use.)"""
        if self.aniPerp0 is None:
            return None
        if d == self.aniShift and self.aniPerp is not None:
            return self.aniPerp
        return self.shiftRows(self.aniPerp0, self.delay, d)

    def aniShiftProposal(self):
        """Delay range proposed for the estimate of the time shift: 25
        delays on either side of the delay closest to 0, where the signal
        rises. Without a delay 0 in the data: the first quarter of the
        delays."""
        t = np.asarray(self.delay, float).ravel()
        n = t.size
        if n == 0:
            return [np.nan, np.nan]
        if t[0] <= 0 <= t[-1]:
            j = int(np.argmin(np.abs(t)))
            return [float(t[max(0, j - 25)]), float(t[min(n - 1, j + 25)])]
        return [float(t[0]),
                float(t[min(n, max(2, int(np.floor(n / 4.0 + 0.5)))) - 1])]

    def aniShiftText(self, E):
        """What an estimate of the time shift (aniShiftEstimate) says, and
        the value for the field of the Anisotropy window (NaN: the field
        is left as it is). Returns (msg, d)."""
        from trspec.kinetics import matlab_round
        u = self.tUnitTxt()
        d = np.nan
        if not E.ok:
            return E.msg, d
        if E.atLimit:
            return ('The best time shift lies at the limit of the search '
                    '(%.4g %s, a quarter of the range): the two data sets '
                    'do not come to match within it. Widen the range, or '
                    'check the two files. The field was left as it is.'
                    % (E.maxShift, u)), d
        # the value with the digits its standard deviation leaves (two
        # of the standard deviation); four digits without one
        if np.isfinite(E.sd) and E.sd > 0:
            nd = int(min(max(1 - np.floor(np.log10(E.sd)), 0), 12))
            d = float(matlab_round(E.d * 10.0 ** nd)) / 10.0 ** nd
            if d == 0:
                d = 0.0                 # (not -0)
            dTxt = ('%.*f %s +/- %.2g (from the scatter between the '
                    'wavelengths)' % (nd, d, u, E.sd))
        else:
            d = float('%.4g' % E.d)
            dTxt = '%.4g %s' % (d, u)
        pct = 100.0 * (1.0 - E.c1 / E.c0) if E.c0 > 0 else 0.0
        tail = ('%d wavelength%s, %d delays; delay step there %.4g %s'
                % (E.nWl, '' if E.nWl == 1 else 's', E.nT, E.step, u))
        if pct < 1:
            # nothing to correct
            return ('No time shift is needed: the best value, %s, lowers '
                    'the sum of squares between the two data sets in the '
                    'range by less than 1 %% (%s). The field was set to 0.'
                    % (dTxt, tail)), 0.0
        return ('Estimated time shift: %s. It lowers the sum of squares '
                'between the two data sets in the range by %.0f %% (%s). '
                'The value is in the field: "Apply" uses it.'
                % (dTxt, pct, tail)), d

    @classmethod
    def aniShiftCost(cls, P, S, t, tMask, d, step):
        """How well the perpendicular data shifted by d, S(t - d), follow
        the parallel data P in the delays tMask. For every wavelength (row)
        the sum of the squared residuals of P = a S(t - d) + b, with the
        factor a and the offset b of that wavelength chosen by least
        squares: the two signals differ by a factor that depends on the
        wavelength (1 + 2 r against 1 - r), and each may sit on its own
        background. Returns one value per wavelength.

        Interpolated values carry less noise than measured ones, the less
        the farther they lie from a measured delay; a comparison of P
        with an interpolated S would therefore prefer shifts of half a
        delay step, whatever the signals do. So both data sets are
        interpolated: the comparison is made four times, with both moved
        together by -3/8, -1/8, 1/8 and 3/8 of the delay step, and the
        four results are averaged. The noise then weighs the same for
        every d.

        P, S: complete rows (no empty value) on the delays t."""
        c = 0.0
        for k in range(4):
            dl = (k - 1.5) / 4.0 * step
            A = cls.shiftRows(P, t, dl)[:, tMask]
            B = cls.shiftRows(S, t, dl + d)[:, tMask]
            A = A - A.mean(axis=1, keepdims=True)
            B = B - B.mean(axis=1, keepdims=True)
            sxx = np.sum(B * B, axis=1)
            sxy = np.sum(A * B, axis=1)
            syy = np.sum(A * A, axis=1)
            ck = syy.copy()
            ok = sxx > 0
            ck[ok] = syy[ok] - sxy[ok] ** 2 / sxx[ok]
            c = c + np.maximum(ck, 0.0)
        return c / 4.0

    @classmethod
    def aniShiftEstimate(cls, P, S, t, tRange):
        """Time shift d of the perpendicular data S that makes S(t - d)
        follow the parallel data P best in the delays tRange (see
        aniShiftCost; the sum over the wavelengths is made smallest). d is
        searched within a quarter of the width of the range on either
        side of 0: 41 equally spaced values, then a golden-section search
        between the neighbors of the best one (60 steps), as the time
        shift of "Subtract data".

        Only wavelengths count whose data are complete (no empty value)
        in the range, widened by the largest shift, one delay step and
        four more delays on both sides: every shift is then judged on
        the same points.

        The estimate rests on one assumption: within the range the two
        signals have the same shape at every wavelength. An anisotropy
        that changes while the signal rises breaks it.

        Returns an object with
          ok       False when nothing could be estimated (then: msg)
          d        the shift
          sd       its standard deviation, from the scatter between the
                   wavelengths: every wavelength pulls d its own way, and
                   the spread of these pulls, with the curvature of the
                   sum of squares, gives the standard deviation (NaN
                   when it cannot be given)
          c0, c1   sum of squares without a shift and at d
          nWl      wavelengths that counted
          nT       delays in the range
          step     delay step in the range (the median)
          maxShift the limit of the search
          atLimit  True when d lies at that limit"""
        from types import SimpleNamespace
        P = np.atleast_2d(np.asarray(P, float))
        S = np.atleast_2d(np.asarray(S, float))
        t = np.asarray(t, float).ravel()
        n = t.size
        R = SimpleNamespace(ok=False, msg='', d=0.0, sd=np.nan, c0=np.nan,
                            c1=np.nan, nWl=0, nT=0, step=np.nan,
                            maxShift=0.0, atLimit=False)
        tr = [float(v) for v in np.ravel(tRange)][:2]
        if len(tr) < 2 or not np.all(np.isfinite(tr)) or tr[0] == tr[1]:
            R.msg = 'Give the delay range of the estimate.'
            return R
        t1, t2 = min(tr), max(tr)
        tMask = (t >= t1) & (t <= t2)
        R.nT = int(np.count_nonzero(tMask))
        if R.nT < 5:
            R.msg = ('The range holds fewer than 5 delays: the time shift '
                     'cannot be estimated there.')
            return R
        step = float(np.median(np.diff(t[tMask])))
        maxShift = 0.25 * (t2 - t1)
        R.step = step
        R.maxShift = maxShift
        j = np.flatnonzero((t >= t1 - maxShift - step)
                           & (t <= t2 + maxShift + step))
        j1 = max(int(j[0]) - 4, 0)
        j2 = min(int(j[-1]) + 4, n - 1)
        Pj = P[:, j1:j2 + 1]
        Sj = S[:, j1:j2 + 1]
        # (a wavelength set to zero by a mask holds nothing)
        rows = (np.all(np.isfinite(Pj), axis=1)
                & np.all(np.isfinite(Sj), axis=1)
                & np.any(Sj != 0, axis=1))
        R.nWl = int(np.count_nonzero(rows))
        if R.nWl < 1:
            R.msg = ('No wavelength has complete data around this range: '
                     'the time shift cannot be estimated there.')
            return R
        Pj = Pj[rows, :]
        Sj = Sj[rows, :]
        tj = t[j1:j2 + 1]
        mj = tMask[j1:j2 + 1]

        def rowCost(x):
            return cls.aniShiftCost(Pj, Sj, tj, mj, x, step)

        def cost(x):
            return float(np.sum(rowCost(x)))

        g = np.linspace(-maxShift, maxShift, 41)
        cg = np.array([cost(x) for x in g])
        i = int(np.argmin(cg))
        lo = float(g[max(i - 1, 0)])
        hi = float(g[min(i + 1, g.size - 1)])
        gr = (np.sqrt(5.0) - 1) / 2
        x1 = hi - gr * (hi - lo)
        x2 = lo + gr * (hi - lo)
        f1 = cost(x1)
        f2 = cost(x2)
        for _ in range(60):
            if f1 <= f2:
                hi, x2, f2 = x2, x1, f1
                x1 = hi - gr * (hi - lo)
                f1 = cost(x1)
            else:
                lo, x1, f1 = x1, x2, f2
                x2 = lo + gr * (hi - lo)
                f2 = cost(x2)
        d = (lo + hi) / 2
        R.d = float(d)
        R.c0 = cost(0.0)
        R.atLimit = bool(abs(d) >= 0.975 * maxShift)
        h = maxShift / 20.0
        cm = rowCost(d - h)
        cc = rowCost(d)
        cp = rowCost(d + h)
        R.c1 = float(np.sum(cc))
        curv = float(np.sum(cp - 2 * cc + cm)) / h ** 2
        if curv > 0 and R.nWl >= 3 and not R.atLimit:
            gl = (cp - cm) / (2 * h)
            R.sd = float(np.sqrt(np.sum((gl - gl.mean()) ** 2)
                                 * R.nWl / (R.nWl - 1.0)) / curv)
        R.ok = True
        return R

    @staticmethod
    def aniTailSpectrum(P, S, t, tailT, Gu):
        """G of the tail matching, wavelength by wavelength. For every
        wavelength the factor between the parallel (P) and the
        perpendicular (S) data in the delays tailT[0] ... tailT[1], found
        as in aniCombine (orthogonal regression through the origin), and
        its uncertainty from the scatter of the n points of this range
        about that line:
        sG = sqrt(sum((P - G S)^2) / (n - 1)) / sqrt(sum(S^2)).
        The noise is taken from the tail itself, not from the whole
        trace: a chirp correction interpolates along the delays and
        lowers the point-to-point noise of the early, closely spaced
        delays, while it hardly touches the tail.
        A check of the one G that is in use (Gu), not used in any
        calculation: G is the same along the spectrum when these values
        scatter as their noise lets expect. Returns a Result with
        G, sG   arrays; NaN where a wavelength has fewer than 5 delays in
                the range or no common signal there
        ok      the wavelengths whose G is determined to 10 %
                (sG <= 0.1 |G|); only those are compared with Gu
        nOk     their number;  nAll: wavelengths with a value
        lo, hi  smallest and largest G among them (NaN: none)
        ratio   RMS of (G - Gu) / sG among them: about 1 when the scatter
                is what the noise explains, larger when G depends on the
                wavelength or the anisotropy has not decayed everywhere.
                (With sG from n points, (G - Gu) / sG follows Student's t
                with n - 1 degrees of freedom, whose mean square is
                (n - 1) / (n - 3): every term is divided by that.) NaN
                with fewer than 3 such wavelengths, or without noise in
                the data."""
        from trspec._util import Result
        P = np.asarray(P, float)
        S = np.asarray(S, float)
        t = np.asarray(t, float).ravel()
        M = P.shape[0]
        T = Result(G=np.full(M, np.nan), sG=np.full(M, np.nan),
                   ok=np.zeros(M, bool), nOk=0, nAll=0, lo=np.nan,
                   hi=np.nan, ratio=np.nan)
        tailT = [float(v) for v in np.ravel(tailT)]
        if len(tailT) != 2 or not np.all(np.isfinite(tailT)):
            return T
        m = (t >= min(tailT)) & (t <= max(tailT))
        if np.count_nonzero(m) < 5:
            return T
        nPt = np.zeros(M)
        for i in range(M):
            a = P[i, m]
            b = S[i, m]
            k = np.isfinite(a) & np.isfinite(b)
            n = int(np.count_nonzero(k))
            if n < 5:
                continue
            a = a[k]
            b = b[k]
            saa = float(np.sum(a ** 2))
            sbb = float(np.sum(b ** 2))
            sab = float(np.sum(a * b))
            if not sab > 0:
                continue
            g = (saa - sbb + np.sqrt((saa - sbb) ** 2 + 4 * sab ** 2)) \
                / (2 * sab)
            T.G[i] = g
            T.sG[i] = np.sqrt(float(np.sum((a - g * b) ** 2)) / (n - 1)) \
                / np.sqrt(sbb)
            nPt[i] = n
        fin = np.isfinite(T.G)
        T.nAll = int(np.count_nonzero(fin))
        with np.errstate(invalid='ignore'):
            T.ok = fin & np.isfinite(T.sG) & (T.sG <= 0.1 * np.abs(T.G))
            # (data without noise: nothing to compare with)
            use = T.ok & (T.sG > 1e-9 * np.abs(T.G))
        T.nOk = int(np.count_nonzero(T.ok))
        if T.nOk > 0:
            T.lo = float(np.min(T.G[T.ok]))
            T.hi = float(np.max(T.G[T.ok]))
        if np.count_nonzero(use) >= 3 and np.isfinite(Gu):
            z2 = ((T.G[use] - Gu) / T.sG[use]) ** 2
            T.ratio = float(np.sqrt(np.mean(
                z2 * (nPt[use] - 3) / (nPt[use] - 1))))
        return T

    @staticmethod
    def aniCombine(P, S, t, view, Gin, tailOn, tailT, thr):
        """One anisotropy view from the corrected parallel (P) and
        perpendicular (S) data. Returns (V, G, nBlank).
        G, the factor of the perpendicular data: Gin, or with tailOn the
        least-squares factor between the two data sets, P ~ G * S, over
        all wavelengths in the delays tailT[0] ... tailT[1] (tail
        matching: there the anisotropy is taken to have decayed). Not the
        ratio of the sums: a transient spectrum has both signs and its
        sum can be close to zero. And least squares in both data sets
        (orthogonal regression through the origin, equal noise in P and
        S): with the noise in S left out, G comes out too small where the
        tail is weak.
        view 'par': P;  'perp': G * S;  'magic': (P + 2 G S) / 3;
        'r': (P - G S) / (P + 2 G S), left empty (NaN) where
        |P + 2 G S| < thr x noise of that sum at the wavelength."""
        from trspec.auto import auto_noise
        P = np.asarray(P, float)
        S = np.asarray(S, float)
        t = np.asarray(t, float).ravel()
        G = float(Gin)
        nBlank = 0
        tailT = [float(v) for v in np.ravel(tailT)]
        if tailOn and len(tailT) == 2 and np.all(np.isfinite(tailT)):
            m = (t >= min(tailT)) & (t <= max(tailT))
            a = P[:, m]
            b = S[:, m]
            ok = np.isfinite(a) & np.isfinite(b)
            saa = float(np.sum(a[ok] ** 2))
            sbb = float(np.sum(b[ok] ** 2))
            sab = float(np.sum(a[ok] * b[ok]))
            if ok.any() and sbb > 0:
                if sab > 0:
                    G = (saa - sbb + np.sqrt((saa - sbb) ** 2
                                             + 4 * sab ** 2)) / (2 * sab)
                else:
                    G = sab / sbb          # no common signal in the range
        if view == 'perp':
            V = G * S
        elif view == 'magic':
            V = (P + 2 * G * S) / 3
        elif view == 'r':
            den = P + 2 * G * S
            with np.errstate(invalid='ignore', divide='ignore'):
                V = (P - G * S) / den
                sg = np.asarray(auto_noise(den, t), float).reshape(-1, 1)
                low = np.abs(den) < thr * sg       # (no noise value: kept)
            nBlank = int(np.count_nonzero(low))
            V[low] = np.nan
            V[~np.isfinite(V)] = np.nan
        else:
            V = P.copy()
        return V, G, nBlank

    # ------------------------------------------------------------------
    # trace of the Fit trace window: one wavelength, or a band
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # position of a band along the delays (Band position window)
    # ------------------------------------------------------------------
    BAND_METHODS = ('moment', 'parabola', 'lognormal')
    BAND_METHOD_NAMES = {'moment': 'first moment',
                         'parabola': 'maximum of a parabola through the top',
                         'lognormal': 'maximum of a fitted log-normal band'}

    def bandSetupDefault(self):
        """Setup of the Band position window for the data shown: the
        whole spectral axis as the band, taken as positive or negative
        after the value of largest size in the data; position in cm-1 by
        the first moment; for intensity data on a wavelength axis the
        signal is converted to one per wavenumber; fit of the position
        from delay 0 (the first delay, if there is no delay 0) to the
        last delay with 2 exponentials."""
        x = np.asarray(self.wavelength, float).ravel()
        t = np.asarray(self.delay, float).ravel()
        D = np.asarray(self.deltaA, float)
        sgn = 'pos'
        A = np.abs(D)
        if np.any(np.isfinite(A)):
            if D.ravel()[int(np.nanargmax(A))] < 0:
                sgn = 'neg'
        tFrom = float(t[0])
        if t[0] <= 0 <= t[-1]:
            tFrom = 0.0
        return {'x1': float(np.min(x)), 'x2': float(np.max(x)), 'sign': sgn,
                'posUnit': 'cm-1',
                'density': bool(self.sigKind == 'int'
                                and self.specUnit != 'cm-1'),
                'method': 'moment', 'frac': 0.8, 'minPct': 5.0,
                'tFrom': tFrom, 'tTo': float(t[-1]), 'nExp': 2,
                'tauStart': [], 'tSpec': float(self.selT)
                if self.selT is not None else float(t[0]), 'logT': True}

    def bandSpectra(self, S):
        """The band of the setup S on the axis its position is given
        on. Returns (v, Y, info): v increasing, Y (len(v) x delays) the
        signal with the band positive, converted to a signal per unit of
        v when S asks for it. Wavelengths without data at any delay
        (masked) are left out. info: n (wavelengths), nSkip (left out),
        lo, hi (their range on the spectral axis), err (a text when
        there is no band)."""
        from trspec.bandpos import band_axis
        x = np.asarray(self.wavelength, float).ravel()
        D = np.asarray(self.deltaA, float)
        info = {'n': 0, 'nSkip': 0, 'lo': np.nan, 'hi': np.nan, 'err': ''}
        x1, x2 = float(S['x1']), float(S['x2'])
        if not (np.isfinite(x1) and np.isfinite(x2)):
            info['err'] = 'Enter both ends of the band.'
            return None, None, info
        rows = (x >= min(x1, x2)) & (x <= max(x1, x2))
        keep = rows & np.any(np.isfinite(D), axis=1)
        info['nSkip'] = int(np.count_nonzero(rows & ~keep))
        info['n'] = int(np.count_nonzero(keep))
        if info['n'] < 3:
            info['err'] = ('The band holds %d wavelength%s with data; at '
                           'least 3 are needed.'
                           % (info['n'], '' if info['n'] == 1 else 's'))
            return None, None, info
        xr = x[keep]
        info['lo'] = float(np.min(xr))
        info['hi'] = float(np.max(xr))
        other = S['posUnit'] != self.specUnit
        if other and np.any(xr <= 0):
            info['err'] = ('The spectral axis has values that are not above '
                           '0: it cannot be converted to %s.' % S['posUnit'])
            return None, None, info
        v, jac = band_axis(xr, self.specUnit, S['posUnit'])
        Y = D[keep, :]
        if other and S['density']:
            Y = Y * jac.reshape(-1, 1)
        if S['sign'] == 'neg':
            Y = -Y
        k = np.argsort(v, kind='stable')
        return v[k], Y[k, :], info

    def bandTrack(self, S):
        """Position of the band of the setup S at every delay
        (trspec.bandpos.band_position). Returns (T, info); T is None
        when there is no band (info['err']). T.v: the position axis."""
        from trspec.bandpos import band_position
        v, Y, info = self.bandSpectra(S)
        if info['err']:
            return None, info
        T = band_position(v, Y, S['method'], float(S['frac']),
                          float(S['minPct']) / 100.0)
        T.v = v
        return T, info

    @staticmethod
    def pmText(x, sd, unit):
        """'value +/- standard deviation unit', both with the decimals
        that show two digits of the standard deviation (more than 12
        decimals: six digits of the value, two of the standard
        deviation). Without a usable standard deviation: the value with
        six digits, marked when the standard deviation is infinite (not
        determined by the data)."""
        if np.isfinite(sd) and sd > 0:
            nd = int(max(1 - np.floor(np.log10(sd)), 0))
            if nd > 12:
                return '%.6g +/- %.2g %s' % (x, sd, unit)
            return '%.*f +/- %.*f %s' % (nd, x, nd, sd, unit)
        if np.isinf(sd):
            return '%.6g %s  [not determined by the data]' % (x, unit)
        return '%.6g %s' % (x, unit)

    def bandShiftLines(self, S, T, info, F):
        """Result text of the Band position window: the band, how its
        position was taken, and the fit F of the position (None: no fit
        was made)."""
        u = self.tUnitTxt()
        pu = S['posUnit']
        left = ''
        if info['nSkip']:
            left = ' (%d without data left out)' % info['nSkip']
        L = ['Band: %.6g ... %.6g %s, %d wavelengths%s, taken as %s'
             % (info['lo'], info['hi'], self.xUnitTxt(), info['n'], left,
                'negative' if S['sign'] == 'neg' else 'positive')]
        how = 'Position: %s, in %s' % (self.BAND_METHOD_NAMES[S['method']],
                                       pu)
        if S['method'] == 'parabola':
            how += ' (points from %.3g of the maximum)' % float(S['frac'])
        L.append(how)
        if pu != self.specUnit:
            if S['density']:
                L.append('Signal: converted to a signal per %s (x %s^2 / 1e7)'
                         % (pu, 'wavenumber' if self.specUnit == 'cm-1'
                            else 'wavelength'))
            else:
                L.append('Signal: as it is (not converted to a signal per '
                         '%s)' % pu)
        L.append('Found at %d of %d delays' % (
            int(np.count_nonzero(T.ok)), int(np.size(T.ok))))
        L.append('')
        if F is None:
            return L
        if not F.ok:
            L.append('No fit of the position: %s' % F.msg)
            return L
        tf = np.asarray(self.delay, float).ravel()[np.asarray(F.mask, bool)]
        L.append('Fit of the position, %.6g ... %.6g %s (%d delays):'
                 % (tf[0], tf[-1], u, F.n))
        for i in range(F.tau.size):
            share = ''
            if F.shift != 0:
                share = '  (%.0f %%)' % (100.0 * F.A[i] / F.shift)
            L.append('  tau_%d = %s' % (i + 1, self.pmText(
                F.tau[i], F.tau_std[i], u)))
            L.append('    A_%d = %s%s' % (i + 1, self.pmText(
                F.A[i], F.A_std[i], pu), share))
        if np.any(np.asarray(F.tau_lim, bool)):
            L.append('  A time constant lies at the limit of the search: the '
                     'data do not determine it, and the values below rest '
                     'on it. Use fewer exponentials.')
        L.append('  at long delays: %s'
                 % self.pmText(F.v_inf, F.v_inf_std, pu))
        L.append('  at delay 0 (extrapolated): %s'
                 % self.pmText(F.v0, F.v0_std, pu))
        L.append('  total shift: %.5g %s' % (F.shift, pu))
        if np.all(F.A > 0) or np.all(F.A < 0):
            L.append('  average time, sum(A tau) / sum(A): %.4g %s'
                     % (F.tau_avg, u))
        L.append('  RMS of the residual: %.3g %s' % (F.rms, pu))
        L.append('')
        L.append('C(t) = (position - position at long delays) / total shift')
        return L

    def bandExportTable(self, S, T, F):
        """Table that the Band position window saves: one row per
        delay; position, height, width, area and asymmetry of the band,
        the fitted position, C(t) and its fit (empty where there is
        none; no C(t) when a time constant of the fit lies at the limit
        of the search). Returns (header, columns)."""
        t = np.asarray(self.delay, float).ravel()
        n = t.size
        fit = np.full(n, np.nan)
        c = np.full(n, np.nan)
        cfit = np.full(n, np.nan)
        if F is not None and F.ok:
            m = np.asarray(F.mask, bool)
            fit[m] = np.asarray(F.fit, float)[m]
            if F.shift != 0 and not np.any(np.asarray(F.tau_lim, bool)):
                c = (np.asarray(T.pos, float) - F.v_inf) / F.shift
                cfit[m] = (fit[m] - F.v_inf) / F.shift
        pu = S['posUnit']
        hdr = ['delay_%s' % self.tUnitHdr(), 'position_%s' % pu, 'height',
               'width_fwhm_%s' % pu, 'area', 'asymmetry', 'fit_%s' % pu,
               'C', 'C_fit']
        cols = [t, np.asarray(T.pos, float), np.asarray(T.height, float),
                np.asarray(T.width, float), np.asarray(T.area, float),
                np.asarray(T.asym, float), fit, c, cfit]
        return hdr, cols

    @staticmethod
    def bandValue(X, x, mode):
        """One trace from the rows X (wavelength x delay) of a band with
        the axis values x. 'mean': mean over the wavelengths;
        'integral': integral over the spectral axis (trapezoid). A delay
        at which one of the wavelengths has no value stays empty, so
        that the same wavelengths contribute at every delay."""
        X = np.asarray(X, float)
        x = np.asarray(x, float).ravel()
        y = np.full(X.shape[1], np.nan)
        ok = np.all(np.isfinite(X), axis=0)
        if not ok.any():
            return y
        if mode == 'mean':
            y[ok] = np.mean(X[:, ok], axis=0)
        else:
            k = np.argsort(x, kind='stable')
            Xs = X[k][:, ok]
            dx = np.diff(x[k]).reshape(-1, 1)
            y[ok] = np.sum(dx * (Xs[:-1, :] + Xs[1:, :]), axis=0) / 2
        return y

    def kinTrace(self, w1, isBand, w2=np.nan, mode='integral'):
        """(d, info): the trace the Fit trace window works on. Not a
        band: the data at the wavelength nearest to w1. A band: the
        wavelengths from w1 to w2, as their integral over the spectral
        axis or as their mean (bandValue). Wavelengths without data at
        any delay (masked) are left out. In the view Anisotropy r the
        band values of the parallel and of the perpendicular data are
        formed first and r is calculated from them: r is a ratio, and
        its mean over a band would be ruled by the wavelengths with
        little signal. info['err'] is a text when there is no trace."""
        from trspec.auto import auto_noise
        x = np.asarray(self.wavelength, float).ravel()
        info = {'band': False, 'wl': np.nan, 'iw': -1, 'w2': np.nan,
                'lo': np.nan, 'hi': np.nan, 'n': 0, 'nSkip': 0, 'mode': '',
                'isR': False, 'err': ''}
        if not isBand:
            iw = int(np.argmin(np.abs(x - w1)))
            info.update(wl=float(x[iw]), iw=iw, n=1)
            return np.array(self.deltaA[iw, :], float), info
        mode = 'mean' if mode == 'mean' else 'integral'
        info.update(band=True, mode=mode, wl=float(w1), w2=float(w2))
        d = np.full(np.size(self.delay), np.nan)
        if not (np.isfinite(w1) and np.isfinite(w2)):
            info['err'] = 'Enter both ends of the band.'
            return d, info
        rows = np.where((x >= min(w1, w2)) & (x <= max(w1, w2)))[0]
        isR = bool(self.aniOn and self.aniView == 'r'
                   and self.aniPar is not None and self.aniPerp is not None)
        X = self.aniPar if isR else self.deltaA
        use = rows[np.any(np.isfinite(X[rows, :]), axis=1)]
        info.update(n=int(use.size), nSkip=int(rows.size - use.size),
                    isR=isR)
        if use.size < (2 if mode == 'integral' else 1):
            info['err'] = ('The band holds fewer than 2 wavelengths with '
                           'data.' if mode == 'integral' else
                           'The band holds no wavelength with data.')
            return d, info
        info.update(lo=float(np.min(x[use])), hi=float(np.max(x[use])))
        if isR:
            P = self.bandValue(self.aniPar[use, :], x[use], mode)
            S = self.bandValue(self.aniPerp[use, :], x[use], mode)
            G = float(self.aniGUsed)
            den = P + 2 * G * S
            with np.errstate(invalid='ignore', divide='ignore'):
                d = (P - G * S) / den
                sg = float(np.ravel(auto_noise(den.reshape(1, -1),
                                               self.delay))[0])
                low = np.abs(den) < self.aniThr * sg   # (no noise value: kept)
            d[low] = np.nan
            d[~np.isfinite(d)] = np.nan
        else:
            d = self.bandValue(X[use, :], x[use], mode)
        return d, info

    def kinTraceLabel(self, info, what):
        """Texts of the Fit trace window for the trace that info (of
        kinTrace) describes. what: 'at' (for plot titles), 'ylabel',
        'file' (suggested name without extension), 'export' (title of
        the file dialog), 'param' (line of the parameter file), 'lines'
        (first lines of the result box)."""
        band = bool(info.get('band'))
        ut = self.xUnitTxt()
        name = 'Wavenumber' if self.specUnit == 'cm-1' else 'Wavelength'
        if what == 'at':
            if band:
                return 'band %.2f ... %.2f %s (%s)' % (
                    info['lo'], info['hi'], self.xUnitAx(), info['mode'])
            return '%s = %.2f %s' % (self.xSym(), info['wl'], self.xUnitAx())
        if what == 'ylabel':
            if band and not info.get('isR'):
                return 'Band %s of %s' % (info['mode'], self.sigName())
            return self.sigName()
        if what == 'file':
            if band:
                return 'kinfit_band_%.1f-%.1f%s' % (info['lo'], info['hi'], ut)
            return 'kinfit_%.1f%s' % (info['wl'], ut)
        if what == 'export':
            if band:
                return 'Export kinetic fit, band %.2f ... %.2f %s' % (
                    info['lo'], info['hi'], ut)
            return 'Export kinetic fit @ %.2f %s' % (info['wl'], ut)
        how = ('%s over %d wavelengths' if info.get('mode') == 'integral'
               else '%s of %d wavelengths') % (info.get('mode'), info.get('n'))
        if what == 'param':
            if not band:
                return '%s: %.4f %s' % (name, info['wl'], ut)
            s = 'Band: %.4f ... %.4f %s, %s' % (info['lo'], info['hi'], ut,
                                                how)
            if info.get('nSkip'):
                s += ' (%d without data left out)' % info['nSkip']
            if info.get('isR'):
                s += ', r from the band values of the two polarizations'
            return s
        if not band:
            return ['%-11s: %.2f %s' % (name, info['wl'], ut)]
        L = ['%-11s: %.2f ... %.2f %s' % ('Band', info['lo'], info['hi'], ut),
             '  ' + how]
        if info.get('nSkip'):
            L.append('  (%d without data left out)' % info['nSkip'])
        if info.get('isR'):
            L += ['  r from the band values of the two', '  polarizations']
        return L

    # ------------------------------------------------------------------
    # log of the steps that make corrected data the starting point (crop,
    # conversion, corrections of Auto Analysis): with it a second data
    # set measured on the axes of the loaded data can be brought to the
    # present stage at any time
    # ------------------------------------------------------------------
    def clearCorrLog(self):
        self.corrLog = []
        self.corrLogOK = True
        self.pendingAuto = None

    def pipeState(self):
        """The corrections in force (background, chirp, masks) as a
        dictionary that pipeOn can apply to any data set on the present
        axes."""
        n = self.wavelength.size
        p = {'bgOn': False, 'bgN': 0, 'bg': None, 'chirpOn': False,
             't0': None, 'nanRows': np.zeros(n, bool),
             'zeroRows': np.zeros(n, bool)}
        if self.bgApplied and self.bgSpectrum is not None:
            p['bgOn'] = True
            p['bgN'] = int(self.bgN)
            p['bg'] = np.asarray(self.bgSpectrum, float).ravel().copy()
        if self.chirpApplied and self.chirpParams is not None:
            p['chirpOn'] = True
            p['t0'] = np.asarray(chirp_model(self.chirpParams,
                                             self.wavelength), float).ravel()
        for wl1, wl2, mode in self.maskedRegions:
            row = (self.wavelength >= wl1) & (self.wavelength <= wl2)
            if mode == 'nan':
                p['nanRows'] |= row
                p['zeroRows'] &= ~row
            else:
                p['zeroRows'] |= row
                p['nanRows'] &= ~row
        return p

    @staticmethod
    def autoLogEntry(P):
        """The corrections of an Auto Analysis run (P) in the form
        pipeOn applies."""
        n = np.asarray(P.wlIn).size
        e = {'bgOn': P.bgStatus == 'applied', 'bgN': int(P.nBg),
             'bg': np.asarray(P.bg, float).ravel().copy(),
             'chirpOn': P.chirp['status'] == 'applied', 't0': None,
             'nanRows': np.zeros(n, bool), 'zeroRows': np.zeros(n, bool)}
        if e['chirpOn']:
            e['t0'] = np.asarray(P.chirp['t0'], float).ravel().copy()
        e['nanRows'][np.asarray(P.maskRows, bool).ravel()] = True
        return e

    @staticmethod
    def pipeOn(X, p, t):
        """Apply the corrections p to another data set X measured on the
        same axes (delays t): its OWN average of the first p['bgN'] delays
        as background (the stored spectrum when no number of delays is
        known), the same time-zero curve, the same masked rows."""
        X = np.array(X, float)
        if p['bgOn']:
            if int(p['bgN']) >= 1:
                seg = X[:, :int(p['bgN'])]
                ok = np.isfinite(seg)
                cnt = ok.sum(axis=1)
                with np.errstate(invalid='ignore', divide='ignore'):
                    b = np.where(ok, seg, 0.0).sum(axis=1) / cnt
                b[~np.isfinite(b)] = 0.0
            else:
                b = np.asarray(p['bg'], float).ravel()
            X = X - b.reshape(-1, 1)
        if p['chirpOn']:
            t = np.asarray(t, float).ravel()
            t0 = np.asarray(p['t0'], float).ravel()
            X = AppState.chirpShift(t, X, t0)
        X[np.asarray(p['nanRows'], bool).ravel(), :] = np.nan
        X[np.asarray(p['zeroRows'], bool).ravel(), :] = 0.0
        return X

    def logRebase(self, wlMask, tMask):
        """Record the step that is about to make corrected data the
        starting point: the corrections in force, then those of an Auto
        Analysis run being applied, then the crop."""
        e = {'op': 'crop', 'pipe': self.pipeState(),
             'auto': self.pendingAuto,
             't': np.array(self.delay, float).ravel(),
             'wlMask': np.array(wlMask, bool).ravel(),
             'tMask': np.array(tMask, bool).ravel()}
        self.pendingAuto = None
        self.corrLog.append(e)

    def replayCorrections(self, X):
        """Bring a data set measured on the axes of the loaded data to
        the stage of deltaA_raw: every logged step is applied to it, in
        order. Returns (X, ok)."""
        X = np.asarray(X, float)
        ok = bool(self.corrLogOK and self.originalDeltaA is not None
                  and X.shape == self.originalDeltaA.shape)
        if not ok:
            return X, False
        for e in self.corrLog:
            if e['op'] == 'crop':
                X = self.pipeOn(X, e['pipe'], e['t'])
                if e.get('auto') is not None:
                    X = self.pipeOn(X, e['auto'], e['t'])
                X = X[np.asarray(e['wlMask'], bool).ravel(), :][
                    :, np.asarray(e['tMask'], bool).ravel()]
            elif e['op'] == 'conv':
                X, _ = self.convertSignal(X, e['from'], e['to'], e['unit'])
        return X, True

    def logHasConversion(self):
        return any(e['op'] == 'conv' for e in self.corrLog)

    def updateRefStage(self):
        """The reference at the stage of deltaA_raw (None when it cannot
        be brought there)."""
        self.refStage = None
        if self.refOrig is None:
            return
        X, ok = self.replayCorrections(self.refOrig)
        if ok and self.deltaA_raw is not None and X.shape == self.deltaA_raw.shape:
            self.refStage = X

    def subAllowed(self):
        """Can a reference be subtracted from the present starting data?
        Returns (tf, why): why is the reason when it cannot."""
        if self.aniOn:
            return False, (
                'Perpendicular data are loaded (anisotropy): a reference '
                'cannot be subtracted here. Subtract it from each '
                'polarization first, save the two results with "Export '
                'Corrected Data...", and load those.')
        if self.solvData:
            return False, (
                'A view without the solvent is shown ("View" of the main '
                'window), made with this reference. Choose the view "Data" '
                'first; the subtraction here can be changed after that.')
        if self.subBaked:
            return False, (
                'A reference was subtracted before the data were cropped, '
                'converted or corrected by Auto Analysis: it is part of the '
                'data now. Press "Reset all" in the main window to change '
                'the subtraction.')
        if not self.corrLogOK:
            return False, (
                'This session was saved by an earlier version, after a crop: '
                'the steps made since loading are not recorded, so a '
                'reference cannot be brought to the same state. Press "Reset '
                'all" in the main window, subtract the reference, then '
                'repeat the corrections.')
        if self.logHasConversion():
            return False, (
                'The data were converted (ΔA <-> ΔT/T) after loading: a '
                'reference cannot be subtracted from them. Press "Reset all" '
                'in the main window, subtract the reference, then convert.')
        return True, ''

    def pipelineBase(self):
        """Data the corrections start from: deltaA_raw, minus the
        reference data set if one is subtracted."""
        base = self.deltaA_raw
        if (self.subApplied and self.refStage is not None
                and self.refStage.shape == base.shape):
            base = base - self.subScale * self.shiftRowsKept(
                self.refStage, self.delay, self.subShift)
        return base

    @staticmethod
    def shiftRows(X, t, d):
        """X(t - d) along the delays (the columns), by interp_runs: a
        cubic spline through the points of every row inside the stretches
        of equally spaced delays, linear between the other delays; beyond
        the two ends the first / last value is held. d > 0 moves the
        features to later delays.

        A spline keeps sharp features and oscillations that linear
        interpolation would flatten; where the delay step jumps it would
        overshoot, so it is used only where the delays are equally
        spaced. Empty values (NaN) stay where they are: the result is
        empty where a delay next to t - d is empty."""
        X = np.asarray(X, float)
        t = np.asarray(t, float).ravel()
        n = t.size
        if d == 0 or n < 2:
            return X.copy()
        X2 = X.reshape(-1, n)
        tq = np.minimum(np.maximum(t - d, t[0]), t[-1])
        Y = np.full(X2.shape, np.nan)
        full = np.all(np.isfinite(X2), axis=1)
        if full.any():
            Y[full, :] = interp_runs(t, X2[full, :], tq)
        if not full.all():
            # interval of every query: t[iq] <= tq <= t[iq + 1]
            iq = np.clip(np.searchsorted(t, tq, side='right') - 1, 0, n - 2)
            onL = tq == t[iq]
            onR = tq == t[iq + 1]
            for ii in np.flatnonzero(~full):
                row = X2[ii, :]
                m = np.isfinite(row)
                if np.count_nonzero(m) < 2:
                    continue
                y = interp_runs(t[m], row[m], tq)[0]
                # empty unless both neighbors are measured (on a delay:
                # unless that delay is)
                ok = (m[iq] | onR) & (m[iq + 1] | onL)
                y[~ok] = np.nan
                Y[ii, :] = y
        return Y.reshape(X.shape)

    @staticmethod
    def shiftRowsKept(X, t, d):
        """shiftRows for the corrections of the main window, which are
        rebuilt at every redraw: the result of the last calls is kept
        (cachedCalc)."""
        if d == 0:
            return np.array(X, dtype=float)
        return cachedCalc('shift', lambda a, b, c: AppState.shiftRows(b, a, float(c)),
                          np.asarray(t, float).ravel(), X, d)

    @staticmethod
    def chirpShift(t, X, t0):
        """The chirp correction of the rows of X (wavelength x delay): row
        i at the delays t + t0[i] (trspec.chirp.chirp_shift). The result of
        the last calls is kept (cachedCalc)."""
        return cachedCalc('chirp', chirp_shift, np.asarray(t, float).ravel(),
                          X, np.asarray(t0, float).ravel())

    @classmethod
    def subCost(cls, S, R, t, tMask, d):
        """Sum of squares of S - s*R(t - d) in the delays tMask, with
        the factor s that makes it smallest. Returns (c, s)."""
        Rd = np.atleast_2d(cls.shiftRows(R, t, d))
        a = np.atleast_2d(np.asarray(S, float))[:, tMask]
        b = Rd[:, tMask]
        ok = np.isfinite(a) & np.isfinite(b)
        a = a[ok]
        b = b[ok]
        den = float(np.sum(b ** 2))
        if a.size == 0 or den <= 0:
            return float(np.sum(a ** 2)), 0.0
        s = float(np.sum(a * b)) / den
        return float(np.sum((a - s * b) ** 2)), s

    @classmethod
    def fitSubtraction(cls, S, R, t, tMask, doShift, maxShift):
        """Factor s and time shift d of the reference R that fit the
        data S best (least squares) in the delays tMask: S ~ s * R(t - d).
        doShift False: d = 0. Otherwise d is searched in
        [-maxShift, maxShift]: 41 equally spaced values, then a
        golden-section search between the neighbors of the best one
        (60 steps). Returns (s, d, c0, c1): c0, c1 = sum of squares
        before and after the subtraction."""
        tMask = np.asarray(tMask, bool)
        a = np.atleast_2d(np.asarray(S, float))[:, tMask]
        c0 = float(np.sum(a[np.isfinite(a)] ** 2))
        d = 0.0
        if doShift and maxShift > 0:
            g = np.linspace(-maxShift, maxShift, 41)
            cg = np.array([cls.subCost(S, R, t, tMask, x)[0] for x in g])
            i = int(np.argmin(cg))
            lo = float(g[max(i - 1, 0)])
            hi = float(g[min(i + 1, g.size - 1)])
            gr = (np.sqrt(5.0) - 1) / 2
            x1 = hi - gr * (hi - lo)
            x2 = lo + gr * (hi - lo)
            f1 = cls.subCost(S, R, t, tMask, x1)[0]
            f2 = cls.subCost(S, R, t, tMask, x2)[0]
            for _ in range(60):
                if f1 <= f2:
                    hi, x2, f2 = x2, x1, f1
                    x1 = hi - gr * (hi - lo)
                    f1 = cls.subCost(S, R, t, tMask, x1)[0]
                else:
                    lo, x1, f1 = x1, x2, f2
                    x2 = lo + gr * (hi - lo)
                    f2 = cls.subCost(S, R, t, tMask, x2)[0]
            d = (lo + hi) / 2
        c1, s = cls.subCost(S, R, t, tMask, d)
        return s, d, c0, c1

    @staticmethod
    def sameAxis(a, b):
        """Two axes hold the same values (to 1e-6 of their span)."""
        a = np.asarray(a, float).ravel()
        b = np.asarray(b, float).ravel()
        if a.size != b.size:
            return False
        if a.size == 0:
            return True
        tol = 1e-6 * max(np.max(np.abs(b)) - np.min(np.abs(b)),
                         np.max(np.abs(b)))
        return bool(np.all(np.abs(a - b) <= tol))

    def setReference(self, wl, t, A, desc):
        """Take a parsed file as the reference data set. Its axes must
        be those of the loaded data, as parsed or with the two axes
        interchanged (as "Transpose axes" does). Returns (ok, msg)."""
        if self.solvData:
            return False, (
                'A view without the solvent is shown ("View" of the main '
                'window), made with the file loaded now. Choose the view '
                '"Data" first, then load another file.')
        wl0 = self.originalWavelength
        t0 = self.originalDelay
        wl = np.asarray(wl, float).ravel()
        t = np.asarray(t, float).ravel()
        A = np.asarray(A, float)
        if not (self.sameAxis(wl, wl0) and self.sameAxis(t, t0)):
            wlT, tT, AT = t.copy(), wl.copy(), A.T.copy()
            iw = np.argsort(wlT, kind='stable')
            it = np.argsort(tT, kind='stable')
            wlT, tT, AT = wlT[iw], tT[it], AT[iw, :][:, it]
            if self.sameAxis(wlT, wl0) and self.sameAxis(tT, t0):
                A = AT
            else:
                return False, (
                    'The axes of this file differ from those of the loaded '
                    'data (%d x %d points here, %d x %d in the data). A '
                    'reference must be measured at the same wavelengths and '
                    'delays.' % (wl.size, t.size, wl0.size, t0.size))
        self.refOrig = A.copy()
        self.refDesc = desc
        self.solvOn = False       # the factor found belongs to the old file
        self.clearSolvDef()
        self.updateRefStage()
        return True, ''

    # ------------------------------------------------------------------
    # solvent coherence (Coherence window): the reference data set is
    # the solvent measured alone; its vibrational coherence is taken out
    # of the residual before the FFT (trspec.coherence: solvent_vc,
    # solvent_shifted, fit_solvent_coherence)
    # ------------------------------------------------------------------
    def solvAllowed(self):
        """Can the solvent coherence be subtracted in the Coherence
        window? Returns (tf, why): why is the reason when it cannot."""
        if self.aniOn:
            return False, ('Perpendicular data are loaded (anisotropy): the '
                           'solvent coherence cannot be subtracted here.')
        if not self.corrLogOK:
            return False, (
                'This session was saved by an earlier version, after a crop: '
                'the steps made since loading are not recorded, so the '
                'solvent data cannot be brought to the same state. Press '
                '"Reset all" in the main window and repeat the corrections.')
        return True, ''

    def solventData(self):
        """The reference (solvent) data set with the corrections in
        force: its own background, the same time-zero curve, the same
        masked rows."""
        return self.pipeOn(self.refStage, self.pipeState(), self.delay)

    def solventInData(self, view=None):
        """What recompute takes out of the data in a solvent view of the
        main window (view: 'osc' | 'all'; without it the view shown,
        solvView):
          'osc'  solvDataA x solvent coherence(t - solvDataShift), made as
                 the Coherence window makes it (solvent_vc, solvent_shifted)
                 in the time window solvDataWin: the oscillating part of
                 the solvent data. Zero outside that window.
          'all'  solvDataA x solvent data(t - solvDataShift) at all delays:
                 the whole solvent signal, with the corrections of the
                 data (solventData) and the time shift of "Subtract
                 data..." (shiftRows). Zero where the solvent is not known.
        On the axes of the data. Zero everywhere in the view 'data' and
        when it cannot be made (no solvent data at the stage of the data,
        fewer than 8 delays in the window). Port of solventInData."""
        from trspec.coherence import solvent_shifted, solvent_vc
        if view is None:
            view = self.solvView
        Y = np.zeros((np.size(self.wavelength), np.size(self.delay)))
        if self.refStage is None or self.refStage.shape != Y.shape:
            return Y
        if view == 'all':
            try:
                Ya = self.solvDataA * self.shiftRowsKept(
                    self.solventData(), self.delay, self.solvDataShift)
                Ya = np.where(np.isfinite(Ya), Ya, 0.0)
                return Ya if Ya.shape == Y.shape else Y
            except Exception:           # noqa: BLE001 (as the MATLAB code)
                return Y
        win = [float(x) for x in np.ravel(self.solvDataWin)]
        if view != 'osc' or len(win) != 2 or not np.all(np.isfinite(win)):
            return Y
        try:
            Vs = solvent_vc(self.solventData(), self.delay, win[0], win[1],
                            self.solvDataPad)
            return self.solvDataA * solvent_shifted(Vs, self.delay,
                                                    self.solvDataShift)
        except Exception:               # noqa: BLE001 (as the MATLAB code)
            return Y

    SOLV_VIEWS = ('data', 'osc', 'all')
    SOLV_VIEW_NAMES = {'data': 'Data', 'osc': 'Solvent coherence subtracted',
                       'all': 'Solvent signal subtracted'}

    def solvViewsOK(self):
        """Are the solvent views of the main window available? They are
        while the solvent coherence is subtracted in the Coherence window
        ("Apply" of its dialog) and that analysis has been made with it:
        its factor, time shift and time window are what the views use.
        Port of solvViewsOK."""
        win = np.ravel(np.asarray(self.solvDataWin, float))
        return bool(self.solvOn and not self.aniOn
                    and self.refStage is not None
                    and self.refStage.shape == (np.size(self.wavelength),
                                                np.size(self.delay))
                    and win.size == 2 and np.all(np.isfinite(win)))

    def setSolvDef(self, A, shift, win, pad):
        """What the solvent views take out of the data: the factor A and
        the time shift the Coherence window subtracts with, its time
        window win = [tMin, tMax] and pad, the largest shift the coherence
        is continued for. Called by that window (view "Data") whenever it
        has made the subtraction. With other values, what is kept for the
        two solvent views describes other data and is dropped. Port of
        setSolvDef."""
        if self.solvData:
            return
        pad = float(np.nanmax([abs(float(pad)), abs(float(shift))]))
        new = [float(A), float(shift), float(min(win)), float(max(win)), pad]
        old = [float(self.solvDataA), float(self.solvDataShift)] \
            + [float(x) for x in np.ravel(self.solvDataWin)] \
            + [float(self.solvDataPad)]
        if len(old) != 5 or not np.array_equal(new, old, equal_nan=True):
            self.solvDataA, self.solvDataShift = new[0], new[1]
            self.solvDataWin = [new[2], new[3]]
            self.solvDataPad = new[4]
            self.solvStore = {}
        self.updateSignalLabel()

    def clearSolvDef(self):
        """No solvent views (nothing is subtracted in the Coherence window
        any more, or its factor does not hold for the data any more). Port
        of clearSolvDef."""
        if self.solvData:
            return
        self.solvDataWin = [np.nan, np.nan]
        self.solvDataPad = 0.0
        self.solvStore = {}
        self.updateSignalLabel()

    def updateSignalLabel(self):
        """(the main window shows the state; nothing to show here)"""

    def getChirpBaseData(self):
        base = self.pipelineBase()
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
        # Anisotropy r: values beyond +-0.5 only occur where the sum of
        # the two signals crosses zero; they would take the contrast
        # from the rest of the map.
        if self.aniOn and self.aniView == 'r':
            cl = min(cl, 0.5)
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

    # ------------------------------------------------------------------
    # signal kind: dA <-> dT/T (transient absorption)
    # ------------------------------------------------------------------
    SIG_KINDS = ('dA', 'dTT', 'int')

    @staticmethod
    def sigNameOf(kind):
        """Name of a signal kind for axis labels: 'dA' | 'dTT' | 'int'
        (intensity: TR-PL, Raman)."""
        return {'dTT': 'ΔT/T', 'int': 'Intensity'}.get(kind, 'ΔA')

    def sigName(self):
        """Name of the signal the working data are."""
        s = self.sigNameOf(self.sigKind)
        if self.aniOn:
            s = {'perp': s + ' (perpendicular)', 'magic': s + ' (magic angle)',
                 'r': 'Anisotropy r'}.get(self.aniView, s + ' (parallel)')
        return s

    @classmethod
    def sigShort(cls, kind):
        """'dA' | 'dTT' | 'int' -> text for controls."""
        return cls.sigNameOf(kind)

    def signalText(self):
        """Text beside the conversion button of the main window."""
        if self.wavelength.size == 0:
            return ''
        if self.aniOn:
            txt = 'Signal: %s   |   two polarizations, G = %.4g%s' % (
                self.sigShort(self.sigKind), self.aniGUsed,
                ' (tail matching)' if self.aniTail else '')
            if self.aniShift != 0:
                txt += ', time shift %.4g %s' % (self.aniShift,
                                                 self.tUnitTxt())
            return txt
        if self.sigKind != self.sigOrigKind:
            txt = 'Signal: %s  (converted from %s)' % (
                self.sigShort(self.sigKind), self.sigShort(self.sigOrigKind))
        else:
            txt = 'Signal: %s' % self.sigShort(self.sigKind)
        # what is subtracted from the data
        if self.subApplied:
            txt += '   |   reference subtracted (scale %.4g)' % self.subScale
        if self.solvData:
            if self.solvView == 'all':
                txt += ('   |   solvent signal subtracted (A = %.4g, all '
                        'delays)' % self.solvDataA)
            else:
                txt += ('   |   solvent coherence subtracted (A = %.4g, '
                        '%.4g ... %.4g %s)'
                        % (self.solvDataA, self.solvDataWin[0],
                           self.solvDataWin[1], self.tUnitTxt()))
        return txt

    @staticmethod
    def convertSignal(A, fromKind, toKind, aUnit):
        """dA <-> dT/T, exact:
            dT/T = 10^(-dA) - 1        dA = -log10(1 + dT/T)
        with dA in OD; aUnit = 'mOD' means the dA values are in units of
        1e-3 OD. dT/T is a fraction (not percent). expm1 / log1p keep
        small signals accurate. dT/T <= -1 has no dA: those values become
        NaN and are counted. Returns (B, nBad)."""
        A = np.asarray(A, float)
        u = 1e-3 if aUnit == 'mOD' else 1.0
        if fromKind == toKind:
            return A.copy(), 0
        if toKind == 'dTT':
            return np.expm1(-np.log(10.0) * u * A), 0
        x = A.copy()
        with np.errstate(invalid='ignore'):
            bad = x <= -1                     # NaN stays NaN below
        nBad = int(np.count_nonzero(bad))
        x[bad] = np.nan
        with np.errstate(invalid='ignore', divide='ignore'):
            B = -np.log1p(x) / (np.log(10.0) * u)
        return B, nBad

    def mapTitle(self, scaleStr=None):
        parts = []
        if self.subApplied:
            parts.append('ref. subtracted')
        if self.bgApplied:
            parts.append('BG(N=%d)' % self.bgN)
        if self.chirpApplied:
            parts.append('chirp(Sellmeier)')
        if self.solvData:
            parts.append('solvent signal subtracted'
                         if self.solvView == 'all'
                         else 'solvent coh. subtracted')
        base = ', '.join(parts) if parts else 'raw'
        if self.isCropped():
            base += ', cropped'
        if self.sigKind != self.sigOrigKind:
            base += ', converted'
        if scaleStr:
            return '2D %s  (%s, %s scale)' % (self.sigName(), base, scaleStr)
        return '2D %s  (%s)' % (self.sigName(), base)

    def statusText(self):
        if self.wavelength.size == 0:
            return 'No data loaded'
        base = self.dataSourceDesc or 'Data'
        s = '%s   (%d λ x %d t)' % (base, self.wavelength.size,
                                        self.delay.size)
        if self.isCropped():
            s += '   [cropped]'
        return s
