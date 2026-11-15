"""
Phase states for comparing maps: a composition set of an ordered phase whose ordered sublattices hold the same site
fractions is the disordered phase. Its label is then the disordered phase's name, or, when that is a disordered part
filtered out of the system (as A1_FCC in the Co-V and Au-Cu databases), the name of a phase of the system with the
same Gibbs energy (FCC_A1 there, which those databases define as an exact copy of the disordered part).

    labeler = Labeler(dbf, comps, phases)
    labeler.label('FCC_4SL', variables, y)       # 'FCC_A1' if disordered, else 'FCC_4SL'
"""
import numpy as np

from pycalphad import calculate

DISORDER_TOL = 1e-4


class Labeler:
    def __init__(self, dbf, comps, phases):
        self.dbf = dbf
        self.ordered = {}       # ordered phase -> (number of ordered sublattices, label of its disordered state)
        for name in phases:
            hints = dbf.phases[name].model_hints
            dis = hints.get('disordered_phase')
            if dis is None or dis == name or hints.get('ordered_phase') != name or dis not in dbf.phases:
                continue
            n_ordered = len(dbf.phases[name].sublattices) - len(dbf.phases[dis].sublattices) + 1
            self.ordered[name] = (n_ordered, self._representative(dis, comps, phases))

    def _representative(self, dis, comps, phases):
        if dis in phases:
            return dis
        structure = [sorted(map(str, c)) for c in self.dbf.phases[dis].constituents]
        candidates = [p for p in phases if p != dis and [sorted(map(str, c)) for c in self.dbf.phases[p].constituents]
                      == structure and list(self.dbf.phases[p].sublattices) == list(self.dbf.phases[dis].sublattices)]
        for p in candidates:
            try:
                ref = calculate(self.dbf, comps, dis, T=[600, 1400], P=101325, N=1, pdens=20, output='GM')
                pts = np.asarray(ref.Y.values).reshape(-1, ref.Y.shape[-1])[:40]
                ref = calculate(self.dbf, comps, dis, T=[600, 1400], P=101325, N=1, points=pts, output='GM')
                other = calculate(self.dbf, comps, p, T=[600, 1400], P=101325, N=1, points=pts, output='GM')
            except Exception:
                continue
            if np.allclose(np.asarray(ref.GM.values), np.asarray(other.GM.values), atol=1e-6, equal_nan=True):
                return p
        return dis

    def disordered(self, phase, variables, y):
        n_ordered = self.ordered[phase][0]
        by_species = {}
        for var, val in zip(variables, y):
            if var.sublattice_index < n_ordered:
                by_species.setdefault(str(var.species), []).append(val)
        return all(max(vals) - min(vals) < DISORDER_TOL for vals in by_species.values())

    def label(self, phase, variables, y):
        if phase in self.ordered and self.disordered(phase, variables, y):
            return self.ordered[phase][1]
        return phase

    def family(self, phase):
        """The phase's disordered state's label for an ordered phase, else its name: the states of one family differ
        only by ordering, which a map does not show where it is of second order (no two-phase field)"""
        return self.ordered[phase][1] if phase in self.ordered else phase

    def label_compset(self, cs):
        pr = cs.phase_record
        return self.label(pr.phase_name, pr.variables, np.asarray(cs.dof)[pr.num_statevars:])
