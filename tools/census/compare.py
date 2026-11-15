"""
Compares mapping census results (census.py outputs).

    python compare.py label1=a.json label2=b.json [label3=c.json ...] [--json summary.json]

Per map and census: ZPF lines, points, failed lines, duplicate pairs (two lines with the same phases where more than
half of the points of the shorter one lie on the other, within 0.05 axis steps and not beyond its ends), zero-width lines (every point's composition sets within 1e-4 of each other),
edge lines (three or more points, all on one edge of the diagram: every composition set lacking the same component),
the phase pairs of the lines and the phases of the nodes. Then, for each pair of censuses, the points of one that
are not on a ZPF line of the other (with the same phases, within half an axis step, in units of the axis steps):
regions one has and the other lacks.
"""
import argparse
import itertools
import json
from collections import Counter

import numpy as np

DUPLICATE_TOL = 0.05      # duplicate pairs: points within this many axis steps of the other line
BEYOND_TOL = 0.01         # ... and not more than this many axis steps beyond its ends
ZERO_WIDTH_TOL = 1e-4
EDGE_TOL = 2e-5           # 2 * MIN_COMPOSITION
COVER_TOL = 0.5           # missing regions: in axis steps


def label(cs):
    """The state of a composition set (states.py): an ordered phase in its disordered state is the disordered phase"""
    return cs.get("label", cs["phase"])


def point_key(point):
    return tuple(sorted(label(cs) for cs in point["sets"]))


def location(point, mp, relative_T):
    """Compositions of all composition sets (by state), and T (relative to the T range, or in T steps)"""
    row = []
    for cs in sorted(point["sets"], key=lambda cs: (label(cs), tuple(cs["X"]))):
        row.extend(np.asarray(cs["X"]) / (1.0 if relative_T else mp["axis"]["x_step"]))
    if mp["axis"]["T"] is not None:
        lo, hi = mp["axis"]["T"]["lims"]
        row.append((point["T"] - lo) / (hi - lo) if relative_T else point["T"] / mp["axis"]["T"]["step"])
    return np.array(row)


def line_runs(line, mp, relative_T):
    """The line's points with all of its phases stable, as runs of consecutive points of the same states:
    [(states, locations)]"""
    runs = []
    for p in line["points"]:
        if len(p["sets"]) != len(line["phases"]):
            continue
        key = point_key(p)
        if not runs or runs[-1][0] != key:
            runs.append((key, []))
        runs[-1][1].append(location(p, mp, relative_T))
    return [(key, np.array(rows)) for key, rows in runs]


def line_locations(line, mp, relative_T):
    rows = [location(p, mp, relative_T) for p in line["points"] if len(p["sets"]) == len(line["phases"])]
    return np.array(rows) if rows else np.zeros((0, 0))


def on_polyline(p, L):
    """Whether p lies on the polyline L (in axis steps) within DUPLICATE_TOL, and not beyond either of its ends"""
    if len(L) == 1:
        return np.max(np.abs(p - L[0])) < DUPLICATE_TOL
    for i, (a, b) in enumerate(zip(L[:-1], L[1:])):
        ab = b - a
        n2 = np.dot(ab, ab)
        t = 0.0 if n2 == 0 else np.dot(p - a, ab) / n2
        if np.linalg.norm(p - (a + np.clip(t, 0, 1) * ab)) >= DUPLICATE_TOL:
            continue
        beyond = (-t if i == 0 and t < 0 else 0.0) + (t - 1 if i == len(L) - 2 and t > 1 else 0.0)
        if beyond * np.sqrt(n2) < BEYOND_TOL:
            return True
    return False


def duplicate_pairs(mp):
    """Pairs of lines with runs of the same states where more than half of the points of the shorter run lie on the
    other (within DUPLICATE_TOL axis steps, and not beyond its ends): the same boundary drawn twice, possibly under
    two names of one state (an ordered phase in its disordered state and the disordered phase)"""
    pairs = set()
    runs = [(i, key, L) for i, z in enumerate(mp["lines"]) for key, L in line_runs(z, mp, False)]
    for (i, ka, La), (j, kb, Lb) in itertools.combinations(runs, 2):
        if i == j or ka != kb or min(len(La), len(Lb)) < 2:
            continue
        A, B = (La, Lb) if len(La) <= len(Lb) else (Lb, La)
        if np.mean([on_polyline(row, B) for row in A]) > 0.5:
            pairs.add((i, j))
    return sorted(pairs)


def zero_width(line):
    pts = [p for p in line["points"] if len(p["sets"]) >= 2]
    if len(pts) < 2:
        return False
    for p in pts:
        X = np.array([cs["X"] for cs in p["sets"]])
        if np.max(X.max(axis=0) - X.min(axis=0)) >= ZERO_WIDTH_TOL:
            return False
    return True


def edge_line(line):
    pts = [p for p in line["points"] if len(p["sets"]) >= 2]
    if len(pts) < 3:
        return False
    for p in pts:
        X = np.array([cs["X"] for cs in p["sets"]])
        if not np.any(np.all(X < EDGE_TOL, axis=0)):
            return False
    return True


def distance_to_polyline(p, L):
    if len(L) == 1:
        return np.linalg.norm(p - L[0])
    a, b = L[:-1], L[1:]
    d = b - a
    dd = np.einsum('ij,ij->i', d, d)
    t = np.clip(np.einsum('ij,ij->i', p - a, d) / np.where(dd > 0, dd, 1), 0, 1)
    return np.min(np.linalg.norm(p - (a + t[:, None] * d), axis=1))


def uncovered(mp_a, mp_b, points=False):
    """Points of mp_a's lines that are not on a line of mp_b with the same states: Counter by states (with
    points=True, also the list of (states, point))"""
    by_states = {}
    for z in mp_b["lines"]:
        for key, L in line_runs(z, mp_b, False):
            by_states.setdefault(key, []).append(L)
    missing = Counter()
    found = []
    for z in mp_a["lines"]:
        for P in z["points"]:
            if len(P["sets"]) != len(z["phases"]):
                continue
            key, p = point_key(P), location(P, mp_a, False)
            others = by_states.get(key, [])
            if not others or min(distance_to_polyline(p, L) for L in others) > COVER_TOL:
                missing[" + ".join(key)] += 1
                found.append((key, P))
    return (missing, found) if points else missing


def summary(mp):
    return {
        "lines": len(mp["lines"]),
        "points": sum(len(z["points"]) for z in mp["lines"]),
        "failed_lines": sum(z["status"] == "FAILED" for z in mp["lines"]),
        "duplicate_pairs": len(duplicate_pairs(mp)),
        "zero_width_lines": sum(zero_width(z) for z in mp["lines"]),
        "edge_lines": sum(edge_line(z) for z in mp["lines"]),
        "phase_pairs": sorted({" + ".join(key) for z in mp["lines"] for key, _ in line_runs(z, mp, False)}),
        "node_phases": sorted({" + ".join(point_key(n)) for n in mp["nodes"]}),
        "elapsed": round(mp["elapsed"], 1),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("censuses", nargs="+", help="label=file.json")
    parser.add_argument("--json", default=None)
    args = parser.parse_args()
    data = {}
    for item in args.censuses:
        label, path = item.split("=", 1)
        with open(path) as f:
            data[label] = {mp["name"]: mp for mp in json.load(f)["maps"]}
    labels = list(data)
    names = [n for n in data[labels[0]] if all(n in data[l] for l in labels)]
    out = {"maps": {}, "uncovered": {}}
    keys = ["lines", "points", "failed_lines", "duplicate_pairs", "zero_width_lines", "edge_lines", "elapsed"]
    print("map".ljust(28) + "".join(f"{l:>34}" for l in labels))
    print(" " * 28 + "".join(f"{'lines/pts/fail/dup/0w/edge/s':>34}" for _ in labels))
    totals = {l: Counter() for l in labels}
    for name in names:
        sums = {l: summary(data[l][name]) for l in labels}
        out["maps"][name] = sums
        row = name[:27].ljust(28)
        for l in labels:
            s = sums[l]
            row += f"{'/'.join(str(s[k]) for k in keys):>34}"
            totals[l].update({k: s[k] for k in keys})
        print(row)
    print("TOTAL".ljust(28) + "".join(f"{'/'.join(str(round(totals[l][k], 1)) for k in keys):>34}" for l in labels))
    print()
    for name in names:
        sums = out["maps"][name]
        base = labels[0]
        for l in labels[1:]:
            for key in ("phase_pairs", "node_phases"):
                gained = sorted(set(sums[l][key]) - set(sums[base][key]))
                lost = sorted(set(sums[base][key]) - set(sums[l][key]))
                if gained or lost:
                    print(f"{name}: {key} {l} vs {base}: +{gained} -{lost}")
    print()
    for a, b in itertools.permutations(labels, 2):
        for name in names:
            missing = uncovered(data[a][name], data[b][name])
            out["uncovered"].setdefault(f"{a} not in {b}", {})[name] = dict(missing)
            if missing:
                print(f"points of {a} not on a line of {b}, {name}: {dict(missing)}")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
