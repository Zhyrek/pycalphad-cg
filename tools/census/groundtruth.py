"""
Ground truth for mapping censuses: where do the maps claim the wrong phases?

    python groundtruth.py plain=plain.json cg=cg.json [--only Co-V,...] [--dbdir DIR] [--json out.json]

For every map, on a grid of probe points across the diagram (binaries: 40 temperatures x 100 compositions, between
the mapping's steps and off simple fractions, where stoichiometric phases are; ternary sections: a 0.02 composition
grid) and at the middle of every tie-line one census draws and another lacks:

1. Each census's claim: the phases its map shows there (inside a three-phase triangle, a two-phase strip between
   consecutive tie-lines, or else the single phase of the nearest tie-line end).
2. The equilibrium with only the claimed phases (the lower Gibbs energy of pycalphad and pycalphad-cg), for every
   claim, and the unrestricted point equilibrium (pycalphad with pdens 1000, and pycalphad-cg).
3. Truth: the lowest of these Gibbs energies (a point equilibrium above it by more than --tol is reported as a miss
   of that solver). A claim is wrong if its Gibbs energy is higher than the truth by more
   than --tol J/mol, or if the equilibrium of its phases has other phases (with phase fractions above 0.05) than it
   claims (a single-phase field claimed where the phase splits, two phases claimed where one is stable). Phases are
   compared by family: an ordered phase and its disordered state are one (a map does not show ordering of second
   order), and a phase that is a copy of the disordered state belongs to it (FCC_A1 beside FCC_4SL).

Errors are grouped by map, census, claimed and true states (as in compare.py: an ordered phase in its disordered state
is the disordered phase).
"""
import argparse
import itertools
import json
import os
import sys
import time
import warnings
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import census  # noqa: E402
from compare import label, uncovered  # noqa: E402
from states import Labeler  # noqa: E402

import pycalphad  # noqa: E402
import pycalphad_cg  # noqa: E402
from pycalphad import variables as v  # noqa: E402
from pycalphad.core.utils import filter_phases, unpack_species  # noqa: E402

NP_MIN = 0.05       # states with smaller phase fractions do not count in structure comparisons
NP_PRESENT = 1e-6


# ----------------------------------------------------------------------------------------------- what a map claims

def _coords(mp, X, T):
    """Plane coordinates of a composition (and T) in axis steps"""
    step = mp['axis']['x_step']
    if mp['ternary']:
        return np.array([X[1] / step, X[2] / step])
    return np.array([X[1] / step, T / mp['axis']['T']['step']])


def _in_triangle(p, a, b, c, eps=1e-9):
    d = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
    if abs(d) < eps:
        return False
    l1 = ((b[1] - c[1]) * (p[0] - c[0]) + (c[0] - b[0]) * (p[1] - c[1])) / d
    l2 = ((c[1] - a[1]) * (p[0] - c[0]) + (a[0] - c[0]) * (p[1] - c[1])) / d
    return min(l1, l2, 1 - l1 - l2) >= -1e-9


class Claims:
    """The phases a census's map shows at a point"""

    def __init__(self, mp, family):
        self.mp, self.family = mp, family
        self.triangles, self.strips, edges, edge_claims = [], [], [], []
        for n in mp['nodes']:
            if mp['ternary'] and len(n['sets']) == 3:
                self.triangles.append(([_coords(mp, cs['X'], n['T']) for cs in n['sets']], self._claim(n['sets'])))
        for z in mp['lines']:
            pts = [p for p in z['points'] if len(p['sets']) == len(z['phases']) == 2]
            for p, q in zip(pts[:-1], pts[1:]):
                if [cs['phase'] for cs in p['sets']] != [cs['phase'] for cs in q['sets']]:
                    continue
                a, b = [_coords(mp, cs['X'], p['T']) for cs in p['sets']]
                c, d = [_coords(mp, cs['X'], q['T']) for cs in q['sets']]
                self.strips.append(((a, b, d, c), self._claim(p['sets'])))
                for (e0, e1), cs in zip([(a, c), (b, d)], p['sets']):   # the phase boundaries: each phase's side
                    edges.append((e0, e1))
                    edge_claims.append(self._claim([cs]))
        self.edges = np.array(edges).reshape(-1, 2, 2)
        self.edge_claims = edge_claims
        n = 16 if mp['ternary'] else 8
        self.rays = np.array([[np.cos(a), np.sin(a)] for a in np.arange(n) * 2 * np.pi / n])

    def _claim(self, sets):
        """(phases to restrict an equilibrium to, the claimed states, the claimed families)"""
        return (tuple(sorted({cs['phase'] for cs in sets})), tuple(sorted(label(cs) for cs in sets)),
                tuple(sorted(self.family(cs['phase']) for cs in sets)))

    def at(self, X, T):
        p = _coords(self.mp, X, T)
        for tri, claim in self.triangles:
            if _in_triangle(p, *tri):
                return claim
        for (a, b, d, c), claim in self.strips:
            if _in_triangle(p, a, b, d) or _in_triangle(p, a, d, c):
                return claim
        return self._single_phase(p)

    def _single_phase(self, p):
        """Outside the two- and three-phase fields: the phase of the first phase boundary met by most of the rays
        from the point (nearest first among equal votes; one ray can slip through the apex of a congruent point
        to the inner boundary)"""
        if not len(self.edges):
            return None
        e0, e = self.edges[:, 0], self.edges[:, 1] - self.edges[:, 0]
        votes = {}
        for d in self.rays:
            # p + t d = e0 + s e  (t > 0, 0 <= s <= 1)
            den = d[0] * e[:, 1] - d[1] * e[:, 0]
            ok = np.abs(den) > 1e-12
            w = e0 - p
            with np.errstate(divide='ignore', invalid='ignore'):
                t = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / den
                s = (w[:, 0] * d[1] - w[:, 1] * d[0]) / den
            hit = ok & (t > 0) & (s >= 0) & (s <= 1)
            if hit.any():
                i = np.flatnonzero(hit)[np.argmin(t[hit])]
                n, dist = votes.get(self.edge_claims[i], (0, np.inf))
                votes[self.edge_claims[i]] = (n + 1, min(dist, t[i]))
        if not votes:
            return None
        return max(votes, key=lambda c: (votes[c][0], -votes[c][1]))


# ----------------------------------------------------------------------------------------------- equilibria

class System:
    def __init__(self, mp, dbf):
        self.mp, self.dbf = mp, dbf
        self.elements = mp['elements']
        self.comps = self.elements + ['VA']
        self.phases = filter_phases(dbf, unpack_species(dbf, self.comps))
        self.labeler = Labeler(dbf, self.comps, self.phases)
        # (a phase that is a copy of an ordered phase's disordered state, see states.py, is of that family)
        copies = {rep: rep for _, rep in self.labeler.ordered.values()}
        self.family = lambda phase: copies.get(phase) or self.labeler.family(phase)
        self.pdens = 500 if mp['ternary'] else 1000

    def grid_conditions(self):
        mp = self.mp
        if mp['ternary']:
            xs = np.arange(0.0097, 1.0, 0.02)
            T = mp['lines'][0]['points'][0]['T'] if mp['lines'] and mp['lines'][0]['points'][0]['T'] else None
            T = T if T is not None else self._section_T()
            return {v.T: T, v.P: 101325, v.N: 1, v.X(self.elements[1]): xs, v.X(self.elements[2]): xs}
        lo, hi = mp['axis']['T']['lims']
        step = (hi - lo) / 40
        return {v.T: lo + step * (np.arange(40) + 0.5), v.P: 101325, v.N: 1,
                v.X(self.elements[1]): np.arange(0.0047, 1.0, 0.01)}

    def _section_T(self):
        for name, db, comps, T, _ in census.TERNARIES:
            if name == self.mp['name']:
                return T
        raise ValueError(self.mp['name'])

    def solve(self, phases, conds, which):
        """Equilibria on a grid (or at a point) with the given phases: list of point results, in the order of
        _points(conds)"""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if which == 'pycalphad':
                ws = pycalphad.Workspace(self.dbf, self.comps, list(phases), conds, calc_opts={'pdens': self.pdens})
            else:
                ws = pycalphad_cg.Workspace(self.dbf, self.comps, list(phases), conds)
            eq = ws.eq
            prf = ws.phase_record_factory
        GM = np.asarray(eq.GM).reshape(-1)
        Ph = np.asarray(eq.Phase).reshape(len(GM), -1)
        NP = np.asarray(eq.NP).reshape(len(GM), -1)
        X = np.asarray(eq.X).reshape(len(GM), Ph.shape[1], -1)
        Y = np.asarray(eq.Y).reshape(len(GM), Ph.shape[1], -1)
        out = []
        for i in range(len(GM)):
            sets = []
            for j in range(Ph.shape[1]):
                name = str(Ph[i, j])
                if not name or not np.isfinite(NP[i, j]) or NP[i, j] < NP_PRESENT:
                    continue
                pr = prf[name]
                y = Y[i, j, :len(pr.variables)]
                sets.append({'phase': name, 'label': self.labeler.label(name, pr.variables, y), 'NP': float(NP[i, j]),
                             'family': self.family(name), 'X': [float(x) for x in X[i, j]]})
            out.append({'GM': float(GM[i]) if np.isfinite(GM[i]) else np.inf, 'sets': _merge(sets)})
        return out

    def points(self, conds):
        """(X, T) of each grid point, in the order of solve()"""
        T = np.atleast_1d(conds[v.T])
        axes = [np.atleast_1d(conds[v.X(el)]) for el in self.elements[1:]]
        out = []
        for t in T:
            for xs in itertools.product(*axes):
                out.append((np.array([1 - sum(xs)] + list(xs)), float(t)))
        return out


def _merge(sets):
    """One composition set per state: sets of the same state and composition (an ordered phase in its disordered
    state beside the disordered phase) are one"""
    out = []
    for cs in sets:
        for o in out:
            if o['label'] == cs['label'] and np.max(np.abs(np.subtract(o['X'], cs['X']))) < 1e-3:
                o['NP'] += cs['NP']
                break
        else:
            out.append(dict(cs))
    return out


def best(results):
    """The lowest-energy result of each point over several solves"""
    return [min(rs, key=lambda r: r['GM']) for rs in zip(*results)]


def states(result, np_min=NP_MIN, key='label'):
    return tuple(sorted(cs[key] for cs in result['sets'] if cs['NP'] >= np_min))


# ----------------------------------------------------------------------------------------------- one map

def check_map(name, maps, dbdirs, tol):
    t0 = time.time()
    first = next(iter(maps.values()))
    dbf = census.load(first['database'])
    system = System(first, dbf)
    claims = {lab: Claims(mp, system.family) for lab, mp in maps.items()}

    # probe points: the grid, then the middle of each disputed tie-line (in batches of its T for binaries)
    conds = system.grid_conditions()
    grid_points = system.points(conds)
    feasible = [i for i, (X, T) in enumerate(grid_points) if X[0] > 1e-9]
    probes = [('grid', grid_points[i]) for i in feasible]
    for a, b in itertools.permutations(maps, 2):
        for key, P in uncovered(maps[a], maps[b], points=True)[1]:
            X = np.mean([cs['X'] for cs in P['sets']], axis=0)
            probes.append((f'{a} only: {" + ".join(key)}', (X, P['T'] if P['T'] is not None else conds[v.T])))

    # claims at every probe
    claim_at = {lab: [c.at(X, T) for _, (X, T) in probes] for lab, c in claims.items()}
    phase_sets = sorted({c[0] for cl in claim_at.values() for c in cl if c is not None})

    # equilibria: unrestricted and restricted to each claimed set of phases, on the grid and at the extra probes
    extra = probes[len(feasible):]
    def solve_all(phases):
        res = []
        for which in ('pycalphad', 'cg'):
            grid = system.solve(phases, conds, which)
            pts = [grid[i] for i in feasible]
            for _, (X, T) in extra:
                c = {v.T: T, v.P: 101325, v.N: 1}
                for el, x in zip(system.elements[1:], X[1:]):
                    c[v.X(el)] = float(np.clip(x, 1e-9, 1))
                pts.append(system.solve(phases, c, which)[0])
            res.append(pts)
        return res
    free = solve_all(system.phases)
    restricted = {ps: best(solve_all(ps)) for ps in phase_sets}
    truth = best(free + list(restricted.values()))
    misses = {which: [i for i in range(len(probes)) if f[i]['GM'] - truth[i]['GM'] > tol]
              for which, f in zip(('pycalphad', 'cg'), free)}

    # verdicts
    errors = defaultdict(list)
    for lab, cl in claim_at.items():
        for i, c in enumerate(cl):
            if c is None:
                continue
            phases, claimed, families = c
            r = restricted[phases][i]
            dG = r['GM'] - truth[i]['GM']
            structure = (states(r, key='family') != families and states(r, 0.0, key='family') != families)
            if dG > tol or structure:
                X, T = probes[i][1]
                errors[(lab, claimed, states(truth[i]))].append(
                    {'kind': probes[i][0], 'X': [round(float(x), 4) for x in X], 'T': round(float(T), 2),
                     'dG': round(float(dG), 3), 'got': states(r, 0.0)})
    return {'name': name, 'probes': len(probes), 'grid': len(feasible), 'elapsed': round(time.time() - t0, 1),
            'solver_misses': {which: [{'X': [round(float(x), 4) for x in probes[i][1][0]],
                                       'T': round(float(probes[i][1][1]), 2),
                                       'dG': round(float(f[i]['GM'] - truth[i]['GM']), 3),
                                       'got': states(f[i], 0.0), 'truth': states(truth[i], 0.0)} for i in misses[which]]
                              for which, f in zip(('pycalphad', 'cg'), free)},
            'errors': [{'census': lab, 'claimed': list(cl), 'truth': list(tr), 'points': pts}
                       for (lab, cl, tr), pts in sorted(errors.items(), key=lambda kv: -len(kv[1]))]}


def describe(err, ternary):
    pts = err['points']
    X = np.array([p['X'] for p in pts])
    T = np.array([p['T'] for p in pts])
    where = (f"x = {X[:, 1].min():.2f}-{X[:, 1].max():.2f}, {X[:, 2].min():.2f}-{X[:, 2].max():.2f}" if ternary
             else f"x = {X[:, 1].min():.3f}-{X[:, 1].max():.3f}, T = {T.min():.0f}-{T.max():.0f} K")
    dG = [p['dG'] for p in pts]
    kinds = sorted({p['kind'] for p in pts})
    return (f"{err['census']}: claims {' + '.join(err['claimed'])}, truth {' + '.join(err['truth'])}: {len(pts)} "
            f"probe(s), {where}, dG up to {max(dG):.2f} J/mol ({', '.join(kinds)})")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('censuses', nargs='+', help='label=file.json')
    parser.add_argument('--only', default=None)
    parser.add_argument('--dbdir', action='append', default=[])
    parser.add_argument('--tol', type=float, default=1.0, help='J/mol')
    parser.add_argument('--json', default=None)
    args = parser.parse_args()
    census.DB_DIRS.extend(args.dbdir)
    data = {}
    for item in args.censuses:
        lab, path = item.split('=', 1)
        with open(path) as f:
            data[lab] = {mp['name']: mp for mp in json.load(f)['maps']}
    labels = list(data)
    names = [n for n in data[labels[0]] if all(n in data[lab] for lab in labels)]
    if args.only:
        names = [n for n in names if n in args.only.split(',')]
    out = []
    for name in names:
        res = check_map(name, {lab: data[lab][name] for lab in labels}, args.dbdir, args.tol)
        out.append(res)
        ternary = data[labels[0]][name]['ternary']
        print(f"== {name}: {res['probes']} probes, {res['elapsed']} s", flush=True)
        for which, ms in res['solver_misses'].items():
            for m in ms[:5]:
                print(f"   {which} point equilibrium above the truth by {m['dG']} J/mol at X={m['X']} T={m['T']}: "
                      f"{' + '.join(m['got'])}, truth {' + '.join(m['truth'])}")
            if len(ms) > 5:
                print(f"   ... {len(ms)} {which} misses")
        for err in res['errors']:
            print('   ' + describe(err, ternary), flush=True)
        if args.json:
            with open(args.json, 'w') as f:
                json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
