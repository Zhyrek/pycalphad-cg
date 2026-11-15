"""Ground truth for census disputes: at the midpoint of each disputed ZPF point's composition sets, the lower-GM
equilibrium of pycalphad (pdens 1000) and pycalphad-cg. Is the disputed line's phase pair stable there?"""
import json, sys, os, warnings, collections
import numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
from compare import line_locations, distance_to_polyline, COVER_TOL
import pycalphad, pycalphad_cg
from pycalphad import Database, variables as v
from pycalphad.core.utils import filter_phases, unpack_species
from importlib.resources import files
HERE = os.path.dirname(os.path.abspath(__file__))
def load(name):
    p = os.path.join(HERE, 'databases', name)
    return Database(p if os.path.exists(p) else str(files('pycalphad.tests.databases').joinpath(name)))
A = {m['name']: m for m in json.load(open(sys.argv[1]))['maps']}
B = {m['name']: m for m in json.load(open(sys.argv[2]))['maps']}
la, lb = sys.argv[3], sys.argv[4]
def disputed(ma, mb):
    by = {}
    for z in mb['lines']:
        L = line_locations(z, mb, False)
        if len(L): by.setdefault(tuple(z['phases']), []).append(L)
    out = []
    for z in ma['lines']:
        for p, P in zip(line_locations(z, ma, False), z['points']):
            o = by.get(tuple(z['phases']), [])
            if not o or min(distance_to_polyline(p, L) for L in o) > COVER_TOL:
                out.append((tuple(z['phases']), P))
    return out
for name in A:
    for (ma, mb, lab) in [(A[name], B[name], la), (B[name], A[name], lb)]:
        pts = disputed(ma, mb)
        if not pts: continue
        dbf = load(ma['database']); els = ma['elements']; comps = els + ['VA']
        phases = filter_phases(dbf, unpack_species(dbf, comps))
        verdict = collections.Counter()
        for pair, P in pts:
            X = np.mean([s['X'] for s in P['sets'][:2]], axis=0)
            T = P['T'] if P['T'] is not None else ma.get('T')
            conds = {v.T: T, v.P: 101325, v.N: 1}
            for el, x in zip(els[1:], X[1:]): conds[v.X(el)] = float(np.clip(x, 1e-9, 1))
            res = []
            for W, kw in [(pycalphad.Workspace, {'calc_opts': {'pdens': 1000}}), (pycalphad_cg.Workspace, {})]:
                eq = W(dbf, comps, phases, conds, **kw).eq
                ph = sorted({str(p) for p, n in zip(eq.Phase.ravel(), eq.NP.ravel()) if p and n > 1e-6})
                res.append((float(np.squeeze(eq.GM)), ph))
            gm, truth = min(res)
            names = set(n.split('#')[0] for n in pair)
            ok = names <= set(truth) and len(truth) <= 3
            verdict[('real' if ok else 'not stable there') + f' (truth {" + ".join(truth)})'] += 1
        print(f'{name}: {lab} only, {" + ".join(pts[0][0])} etc. ({len(pts)} pts):', dict(verdict))
