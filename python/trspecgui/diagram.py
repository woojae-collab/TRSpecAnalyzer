"""Drawing of a kinetic scheme: where the boxes, arrows and labels go.
Port of targetDiagramLayout and tdBoxExit (indices count from 0 here;
the ground-state sink is -1)."""
import math

import numpy as np


def _mround(x):
    """MATLAB round: halves away from zero."""
    return math.floor(x + 0.5) if x >= 0 else -math.floor(-x + 0.5)


def boxExit(P, Q, R, pad):
    """Point where the way from P (inside the box R = [x, y, w, h],
    enlarged by pad) toward Q leaves the box."""
    x0, x1 = R[0] - pad, R[0] + R[2] + pad
    y0, y1 = R[1] - pad, R[1] + R[3] + pad
    dx, dy = Q[0] - P[0], Q[1] - P[1]
    t = 1.0
    if dx > 0:
        t = min(t, (x1 - P[0]) / dx)
    elif dx < 0:
        t = min(t, (x0 - P[0]) / dx)
    if dy > 0:
        t = min(t, (y1 - P[1]) / dy)
    elif dy < 0:
        t = min(t, (y0 - P[1]) / dy)
    t = max(t, 0.0)
    return [P[0] + t * dx, P[1] + t * dy]


def _text(x, y, s, ha, va, size, kind, sp=-1, bg=False):
    return {'x': x, 'y': y, 'str': s, 'ha': ha, 'va': va, 'size': size,
            'kind': kind, 'sp': sp, 'bg': bg, 'bb': None}


def targetDiagramLayout(sch, W, H, fs0=11.0, tauStep=None, undet=None,
                        unitTxt=''):
    """Layout of the scheme sch for an area of W x H (pixels), letters
    of at most fs0 pixels. tauStep: time constant of every step (shown
    under the rate, with the share of each branch); undet: steps whose
    rate the data do not determine (marked '?').

    Returns a dict: W, H (area needed, at least the one given), fs,
    fits, rect (pos = [x, y, w, h], sp, init), line and head (x, y),
    text (x, y, str, ha, va, size, kind, sp, bg, bb), col, row, noteY,
    noteSize. Compartments in the same connected part are set in
    columns by their place along the scheme; parts stand side by side.
    """
    nS = sch.n_s
    st = sch.steps
    nSt = len(st)
    tauStep = [] if tauStep is None else \
        [float(v) for v in np.asarray(tauStep, float).ravel()]
    hasVal = nSt > 0 and len(tauStep) == nSt
    undet = [] if undet is None else [bool(v) for v in np.ravel(undet)]
    if len(undet) != nSt:
        undet = [False] * nSt
    D = {'W': W, 'H': H, 'fs': fs0, 'fits': True, 'rect': [], 'line': [],
         'head': [], 'text': [], 'col': [], 'row': [], 'noteY': 0.0,
         'noteSize': fs0}
    if nS < 1:
        return D
    fromV = [s['from'] for s in st]
    toV = [s['to'] for s in st]
    c0 = np.asarray(sch.c0, float).ravel()
    sinkName = getattr(sch, 'sink_name', '') or 'GS'
    species = list(sch.species)

    # ---- labels of the steps -------------------------------------------
    with np.errstate(divide='ignore', invalid='ignore'):
        kAll = 1.0 / np.asarray(tauStep, float) if hasVal else None
    lab = []
    for q in range(nSt):
        s1 = st[q]['expr']
        if st[q]['coef'] != 1:
            s1 = '%s [x%g]' % (s1, st[q]['coef'])
        one = [s1]
        if hasVal:
            tv = tauStep[q]
            if not (np.isfinite(tv) and tv > 0):
                s2 = '?'
            elif 1e3 <= tv < 1e5:
                s2 = ('%.0f %s' % (tv, unitTxt)).strip()
            else:
                s2 = ('%.3g %s' % (tv, unitTxt)).strip()
            if undet[q]:
                s2 += ' ?'
            outI = [j for j in range(nSt) if fromV[j] == fromV[q]]
            kOut = kAll[outI]
            if len(outI) > 1 and np.all(np.isfinite(kOut)) and \
                    np.sum(kOut) > 0:
                s2 = '%s (%.0f%%)' % (s2, 100.0 / tauStep[q] / np.sum(kOut))
            one.append(s2)
        lab.append(one)

    # ---- remarks inside the boxes --------------------------------------
    defC0 = np.zeros(nS)
    defC0[0] = 1
    showInit = not np.array_equal(c0, defC0)
    nonneg = np.asarray(sch.nonneg, bool).ravel()
    ann = []
    for q in range(nS):
        it = []
        if showInit and c0[q] > 0:
            it.append('init %g' % c0[q])
        nz = 0
        for z in sch.zero:
            if z['sp'] != q:
                continue
            nz += 1
            if nz == 1:
                it.append('zero %g-%g' % (z['lo'], z['hi']))
            elif nz == 2:
                it[-1] = it[-1] + ' ..'
        for r in sch.rel:
            if r['dep'] != q:
                continue
            if r['ratio'] == 1:
                it.append('= ' + species[r['master']])
            else:
                it.append('= %g %s' % (r['ratio'], species[r['master']]))
        if nonneg[q]:
            it.append('nonneg')
        a = ', '.join(it)
        if len(a) > 18:
            a = a[:16] + '..'
        ann.append(a)

    # ---- connected parts, direction along the scheme, layers -----------
    A = np.zeros((nS, nS), bool)
    for q in range(nSt):
        if toV[q] >= 0:
            A[fromV[q], toV[q]] = True
    U = A | A.T
    comp = [0] * nS
    nC = 0
    for s0 in range(nS):
        if comp[s0] > 0:
            continue
        nC += 1
        comp[s0] = nC
        stk = [s0]
        while stk:
            u = stk.pop()
            nb = [v for v in range(nS) if U[u, v] and comp[v] == 0]
            for v in nb:
                comp[v] = nC
            stk += nb
    Adag = A.copy()                  # steps that close a loop are left out
    state = [0] * nS
    order = [i for i in range(nS) if c0[i] > 0] + \
            [i for i in range(nS) if not (c0[i] > 0)]
    for r in order:
        if state[r] != 0:
            continue
        stk, nxt = [r], [0]
        state[r] = 1
        while stk:
            u = stk[-1]
            nb = [v for v in range(nS) if A[u, v]]
            k = nxt[-1]
            if k >= len(nb):
                state[u] = 2
                stk.pop()
                nxt.pop()
            else:
                nxt[-1] = k + 1
                v = nb[k]
                if state[v] == 1:
                    Adag[u, v] = False
                elif state[v] == 0:
                    state[v] = 1
                    stk.append(v)
                    nxt.append(0)
    layer = [0] * nS
    indeg = [int(Adag[:, v].sum()) for v in range(nS)]
    queue = [v for v in range(nS) if indeg[v] == 0]
    while queue:
        u = queue.pop(0)
        for v in [v for v in range(nS) if Adag[u, v]]:
            layer[v] = max(layer[v], layer[u] + 1)
            indeg[v] -= 1
            if indeg[v] == 0:
                queue.append(v)
    for u in [u for u in range(nS) if int(Adag[:, u].sum()) == 0]:
        nb = [v for v in range(nS) if Adag[u, v]]
        if nb:
            layer[u] = min(layer[v] for v in nb) - 1

    E = []
    for u in range(nS):
        for v in range(u + 1, nS):
            if not U[u, v]:
                continue
            lo, hi = (u, v) if layer[u] <= layer[v] else (v, u)
            E.append({
                'lo': lo, 'hi': hi,
                'fwd': [q for q in range(nSt)
                        if fromV[q] == lo and toV[q] == hi],
                'bwd': [q for q in range(nSt)
                        if fromV[q] == hi and toV[q] == lo],
                'dum': []})
    nE = len(E)
    eComp = [comp[e['lo']] for e in E]

    # ---- rows: order inside every column -------------------------------
    ysp = [0.0] * nS
    colOf = [0] * nS
    nRow = 1
    colBase = 0
    colComp, colN = [], []
    for c in range(1, nC + 1):
        mem = [i for i in range(nS) if comp[i] == c]
        nL = max(layer[i] for i in mem) + 1
        ec = [k for k in range(nE) if eComp[k] == c]
        for ly in range(nL):
            it = []
            for s in [i for i in mem if layer[i] == ly]:
                yp = []
                for k in ec:
                    if E[k]['hi'] != s:
                        continue
                    if layer[E[k]['lo']] == ly - 1:
                        yp.append(ysp[E[k]['lo']])
                    elif E[k]['dum']:
                        yp.append(E[k]['dum'][-1])
                key = -math.inf
                if yp:
                    key = sum(yp) / len(yp)
                it.append((-key, 0, s))
            for k in ec:
                if not (layer[E[k]['lo']] < ly < layer[E[k]['hi']]):
                    continue
                if layer[E[k]['lo']] == ly - 1:
                    key = ysp[E[k]['lo']]
                else:
                    key = E[k]['dum'][-1]
                it.append((-key, 1, k))
            it.sort()
            n = len(it)
            for w, row in enumerate(it):
                yv = (n - 1) / 2.0 - w
                if row[1] == 0:
                    ysp[row[2]] = yv
                else:
                    E[row[2]]['dum'].append(yv)
            nRow = max(nRow, n)
            colN.append(n)
        for i in mem:
            colOf[i] = colBase + layer[i]
        colComp += [c] * nL
        colBase += nL
    nCol = colBase

    sideR = [False] * nS
    sideL = [False] * nS
    for e in E:
        yr = [ysp[e['lo']]] + e['dum'] + [ysp[e['hi']]]
        isP = bool(e['fwd']) and bool(e['bwd'])
        if isP or yr[1] < yr[0]:
            sideR[e['lo']] = True
        if isP or yr[-2] < yr[-1]:
            sideL[e['hi']] = True
    eF = [sum((lab[q] for q in e['fwd']), []) for e in E]
    eB = [sum((lab[q] for q in e['bwd']), []) for e in E]
    sinkOf = [[q for q in range(nSt) if fromV[q] == s and toV[q] < 0]
              for s in range(nS)]
    sLab = [sum((lab[q] for q in sinkOf[s]), []) for s in range(nS)]
    anySink = any(bool(x) for x in sinkOf)
    under = [bool(sinkOf[q]) and sideR[q] and sideL[q] for q in range(nS)]
    anyAnn = any(bool(a) for a in ann)
    maxCh = maxLn = maxSCh = maxSLn = 0
    anyPair = False
    for k in range(nE):
        maxCh = max([maxCh] + [len(s) for s in eF[k] + eB[k]])
        maxLn = max(maxLn, len(eF[k]), len(eB[k]))
        anyPair = anyPair or (bool(eF[k]) and bool(eB[k]))
    for q in range(nS):
        maxSCh = max([maxSCh] + [len(s) for s in sLab[q]])
        maxSLn = max(maxSLn, len(sLab[q]))

    # ---- sizes: letters are made smaller until the drawing fits --------
    fs = float(fs0)
    for _pass in range(6):
        fsS = max(6.0, _mround(1.68 * fs) / 2.0)
        cwN = 0.66 * fs
        cwS = 0.56 * fsS
        lh = 1.22 * fsS
        bh = 1.9 * fs
        if anyAnn:
            bh = bh + 1.15 * fsS
        colW = [3.0 * fs] * nCol
        for q in range(nS):
            wq = max(3.0 * fs, len(species[q]) * cwN + 1.3 * fs,
                     len(ann[q]) * cwS + 1.0 * fs)
            colW[colOf[q]] = max(colW[colOf[q]], wq)
        gx = max(2.6 * fs, max(maxCh, maxSCh) * cwS + 1.3 * fs)
        sOver = 0.0
        if anySink:
            sOver = max(0.0, 0.45 * fs + maxSCh * cwS - min(colW) / 2.0)
        gxC = 1.6 * fs + sOver
        mS = 0.8 * fs + sOver
        stubLen = stubH = 0.0
        if anySink:
            stubLen = max(1.7 * fs, maxSLn * lh + 0.6 * fs)
            stubH = stubLen + 1.3 * fsS
            if any(under):
                stubH = stubH + maxSLn * lh
        labH = 0.0
        if nE > 0:
            labH = maxLn * lh + 0.35 * fs
        vGap = 1.0 * fs
        if nRow > 1:
            vGap = max(vGap, 0.5 * labH + 0.5 * fs)
        ry = bh + stubH + vGap
        dOff = 0.34 * fs
        noteH = lh + 0.5 * fs
        mT = 0.6 * fs + max(0.0, labH + dOff - bh / 2.0)
        mB = 0.4 * fs + noteH + max(
            stubH, (1.0 if anyPair else 0.0) * max(0.0,
                                                   labH + dOff - bh / 2.0))
        xc = [0.0] * nCol
        x = mS
        for j in range(nCol):
            if j > 0:
                x += gx if colComp[j] == colComp[j - 1] else gxC
            xc[j] = x + colW[j] / 2.0
            x += colW[j]
        Wn = x + mS
        Hn = mB + nRow * bh + (nRow - 1) * (ry - bh) + mT
        if (Wn <= W and Hn <= H) or fs <= 6:
            break
        sc = min(W / Wn, H / Hn, 0.97)
        fs = max(6.0, math.floor(2.0 * fs * sc) / 2.0)

    ox = max(0.0, (W - Wn) / 2.0)
    oy = max(0.0, (H - Hn) / 2.0)

    def yOf(r):
        return oy + mB + bh / 2.0 + ((nRow - 1) / 2.0 + r) * ry

    xc = [v + ox for v in xc]
    R = []
    T = []
    for q in range(nS):
        wq = colW[colOf[q]]
        Rq = [xc[colOf[q]] - wq / 2.0, yOf(ysp[q]) - bh / 2.0, wq, bh]
        R.append(Rq)
        D['rect'].append({'pos': list(Rq), 'sp': q, 'init': bool(c0[q] > 0)})
        cx = Rq[0] + wq / 2.0
        if not ann[q]:
            T.append(_text(cx, Rq[1] + bh / 2.0, species[q], 'center',
                           'middle', fs, 'name', q))
        else:
            T.append(_text(cx, Rq[1] + bh - 0.95 * fs, species[q], 'center',
                           'middle', fs, 'name', q))
            T.append(_text(cx, Rq[1] + 0.3 * fs + 0.5 * fsS, ann[q],
                           'center', 'middle', fsS, 'annot', q))

    # ---- arrows between compartments -----------------------------------
    ah = 0.8 * fs
    aw = 0.3 * fs
    for k, e in enumerate(E):
        lo, hi = e['lo'], e['hi']
        P = [[R[lo][0] + R[lo][2] / 2.0, R[lo][1] + bh / 2.0]]
        for w in range(len(e['dum'])):
            P.append([xc[colOf[lo] + 1 + w], yOf(e['dum'][w])])
        P.append([R[hi][0] + R[hi][2] / 2.0, R[hi][1] + bh / 2.0])
        isPair = bool(eF[k]) and bool(eB[k])
        for dirn in (1, 2):
            ln = eF[k] if dirn == 1 else eB[k]
            if not ln:
                continue
            sh = dOff * (3 - 2 * dirn) if isPair else 0.0
            Q = [[p[0], p[1] + sh] for p in P]
            rs, re = R[lo], R[hi]
            if dirn == 2:
                Q = Q[::-1]
                rs, re = R[hi], R[lo]
            Q[0] = boxExit(Q[0], Q[1], rs, 0.0)
            Q[-1] = boxExit(Q[-1], Q[-2], re, 0.12 * fs)
            ux, uy = Q[-1][0] - Q[-2][0], Q[-1][1] - Q[-2][1]
            nu = math.sqrt(ux * ux + uy * uy)
            if nu > 0:
                ux, uy = ux / nu, uy / nu
            else:
                ux, uy = 1.0, 0.0
            tip = list(Q[-1])
            bas = [tip[0] - ah * ux, tip[1] - ah * uy]
            nvx, nvy = -uy, ux
            D['head'].append({
                'x': [tip[0], bas[0] + aw * nvx, bas[0] - aw * nvx],
                'y': [tip[1], bas[1] + aw * nvy, bas[1] - aw * nvy]})
            Q[-1] = [tip[0] - 0.6 * ah * ux, tip[1] - 0.6 * ah * uy]
            D['line'].append({'x': [p[0] for p in Q], 'y': [p[1] for p in Q]})
            nL1 = len(ln)
            hw = max(len(s) for s in ln) * cwS / 2.0
            hh = nL1 * lh / 2.0
            onPatch = False
            side = 1.0
            if isPair and dirn == 2:
                side = -1.0
            if not e['dum']:
                M = [(Q[0][0] + tip[0]) / 2.0, (Q[0][1] + tip[1]) / 2.0]
                dv = [tip[0] - Q[0][0], tip[1] - Q[0][1]]
                if abs(dv[1]) < 0.2 * abs(dv[0]):
                    C = [M[0], M[1] + side * (0.3 * fs + hh)]
                else:
                    if not isPair:
                        yRef = 0.0
                        if colN[colOf[lo]] < colN[colOf[hi]]:
                            yRef = ysp[lo]
                        elif colN[colOf[hi]] < colN[colOf[lo]]:
                            yRef = ysp[hi]
                        if (ysp[lo] + ysp[hi]) / 2.0 < yRef:
                            side = -1.0
                    nd = math.sqrt(dv[0] ** 2 + dv[1] ** 2)
                    nx, ny = -dv[1] / nd, dv[0] / nd
                    if ny < 0:
                        nx, ny = -nx, -ny
                    f = side * (hw * abs(nx) + hh * abs(ny) + 0.3 * fs)
                    C = [M[0] + f * nx, M[1] + f * ny]
            else:
                iv = int(math.ceil(len(e['dum']) / 2.0))
                V = [P[iv][0], P[iv][1] + sh]
                isLow = P[iv][1] < P[iv - 1][1] and P[iv][1] < P[iv + 1][1]
                isHigh = P[iv][1] > P[iv - 1][1] and P[iv][1] > P[iv + 1][1]
                if isPair:
                    C = [V[0], V[1] + side * (0.35 * fs + hh)]
                    onPatch = True
                elif isLow:
                    C = [V[0], V[1] - (0.35 * fs + hh)]
                elif isHigh:
                    C = [V[0], V[1] + (0.35 * fs + hh)]
                else:
                    C = V
                    onPatch = True
            for w in range(nL1):
                knd = 'value' if hasVal and (w + 1) % 2 == 0 else 'rate'
                T.append(_text(C[0], C[1] + ((nL1 - 1) / 2.0 - w) * lh,
                               ln[w], 'center', 'middle', fsS, knd,
                               bg=onPatch))

    # ---- losses to the ground state ------------------------------------
    for q in range(nS):
        if not sinkOf[q]:
            continue
        cx = R[q][0] + R[q][2] / 2.0
        yb = R[q][1]
        tipY = yb - stubLen
        D['line'].append({'x': [cx, cx], 'y': [yb, tipY + 0.6 * ah]})
        D['head'].append({'x': [cx, cx + aw, cx - aw],
                          'y': [tipY, tipY + ah, tipY + ah]})
        T.append(_text(cx, tipY - 0.15 * fs, sinkName, 'center', 'top', fsS,
                       'sink'))
        ln = sLab[q]
        for w in range(len(ln)):
            knd = 'value' if hasVal and (w + 1) % 2 == 0 else 'rate'
            ty = yb - 0.25 * fs - w * lh
            if under[q]:
                tx, tha = cx, 'center'
                ty = tipY - 0.15 * fs - 1.3 * fsS - w * lh
            elif sideR[q]:
                tx, tha = cx - 0.45 * fs, 'right'
            else:
                tx, tha = cx + 0.45 * fs, 'left'
            T.append(_text(tx, ty, ln[w], tha, 'top', fsS, knd))

    note = 'thick frame = initial population'
    if hasVal and any(undet):
        note += ';  ? = not determined'
    D['noteY'] = 0.4 * fs
    D['noteSize'] = fsS
    T.append(_text(0.6 * fs, D['noteY'], note, 'left', 'bottom', fsS, 'note'))
    for t in T:
        cw = 0.66 * t['size'] if t['kind'] == 'name' else 0.56 * t['size']
        tw = len(t['str']) * cw
        th = 1.15 * t['size']
        if t['ha'] == 'center':
            x0 = t['x'] - tw / 2.0
        elif t['ha'] == 'right':
            x0 = t['x'] - tw
        else:
            x0 = t['x']
        if t['va'] == 'middle':
            y0 = t['y'] - th / 2.0
        elif t['va'] == 'top':
            y0 = t['y'] - th
        else:
            y0 = t['y']
        t['bb'] = [x0, y0, x0 + tw, y0 + th]
    D['text'] = T
    D['fs'] = fs
    D['fits'] = bool(Wn <= W + 1e-6 and Hn <= H + 1e-6)
    D['W'] = max(W, Wn)
    D['H'] = max(H, Hn)
    D['col'] = list(colOf)
    D['row'] = list(ysp)
    return D
