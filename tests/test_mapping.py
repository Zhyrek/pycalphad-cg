"""Mapping with column generation through pycalphad's own strategies."""
import matplotlib
import numpy as np
import pycalphad
import pytest

import pycalphad_cg as pyc
from conftest import load_database

matplotlib.use('Agg')


def _conditions(element):
    return {pyc.v.T: (1000, 2000, 25), pyc.v.X(element): (0, 1, 0.025), pyc.v.P: 101325, pyc.v.N: 1}


def test_binary_strategy_uses_column_generation(cg_calls):
    dbf = load_database('alni_dupin_2001.tdb')
    strategy = pyc.mapping.BinaryStrategy(dbf, ['AL', 'NI', 'VA'], list(dbf.phases), _conditions('AL'))
    strategy.do_map()
    assert cg_calls['workspace'] > 0 and cg_calls['local'] > 0 and cg_calls['global_min'] > 0
    assert len(strategy.zpf_lines) > 0
    # the two-phase field next to pure Ni, FCC_A1 (as FCC_L12) + L12 at 1000 K: compare with an equilibrium
    tielines = {tuple(sorted(line.stable_phases)) for line in strategy.zpf_lines}
    assert any('FCC_L12' in phases for phases in tielines)
    calls = dict(cg_calls)
    plain = pycalphad.mapping.BinaryStrategy(dbf, ['AL', 'NI', 'VA'], list(dbf.phases), _conditions('AL'))
    plain.do_map()
    assert cg_calls == calls   # pycalphad's own strategy is unaffected


def test_binplot(cg_calls):
    dbf = load_database('alfe.tdb')
    ax = pyc.binplot(dbf, ['AL', 'FE', 'VA'], list(dbf.phases), _conditions('AL'))
    assert cg_calls['local'] > 0
    assert len(ax.get_lines()) + len(ax.collections) > 0


def test_ternplot(cg_calls):
    dbf = load_database('alcrni.tdb')
    conds = {pyc.v.T: 1273, pyc.v.X('AL'): (0, 1, 0.025), pyc.v.X('CR'): (0, 1, 0.025), pyc.v.P: 101325,
             pyc.v.N: 1}
    ax = pyc.ternplot(dbf, ['AL', 'CR', 'NI', 'VA'], list(dbf.phases), conds)
    assert cg_calls['local'] > 0
    assert len(ax.get_lines()) + len(ax.collections) > 0


def test_step_strategy(cg_calls):
    dbf = load_database('alfe.tdb')
    conds = {pyc.v.T: (600, 1800, 20), pyc.v.X('AL'): 0.3, pyc.v.P: 101325, pyc.v.N: 1}
    strategy = pyc.mapping.StepStrategy(dbf, ['AL', 'FE', 'VA'], list(dbf.phases), conds)
    strategy.do_map()
    assert cg_calls['local'] > 0
    plain = pycalphad.mapping.StepStrategy(dbf, ['AL', 'FE', 'VA'], list(dbf.phases), conds)
    plain.do_map()
    # the same phases appear along the step
    phases = [sorted(s.get_data(pyc.v.T, pyc.v.X('AL')).phases) for s in (strategy, plain)]
    assert phases[0] == phases[1] and len(phases[0]) > 0


def test_step_across_a_closing_miscibility_gap():
    "A two-phase step above the critical point ends in one composition set, as with pycalphad's Solver."
    import copy
    import pycalphad.mapping.zpf_equilibrium as zeq
    from pycalphad.mapping.starting_points import point_from_equilibrium
    dbf = load_database('alzn_mey.tdb')
    conds = {pyc.v.T: 620, pyc.v.P: 101325, pyc.v.N: 1, pyc.v.X('AL'): 0.6}
    with pyc.cg_mode():
        point = point_from_equilibrium(dbf, ['AL', 'ZN', 'VA'], ['HCP_A3', 'LIQUID', 'FCC_A1'], copy.deepcopy(conds))
        assert [cs.phase_record.phase_name for cs in point.stable_composition_sets] == ['FCC_A1', 'FCC_A1']
        result = zeq.update_equilibrium_with_new_conditions(point, {**conds, pyc.v.T: 630})
    assert result is not None
    sets = result[0].stable_composition_sets
    assert len(sets) == 1 and abs(float(sets[0].X[0]) - 0.6) < 1e-9


def test_flat_energy_without_mixing_entropy_is_converged_precisely():
    "Final results are converged in the site fractions also where the energy has little curvature (no ideal mixing)."
    from pycalphad import Database, Model, variables as v
    dbf = Database("""
    ELEMENT A A 1 1 1 !
    ELEMENT B B 1 1 1 !
    TYPE_DEFINITION % SEQ * !
    DEFINE_SYSTEM_DEFAULT VA !
    PHASE S %  1 1.0 !
    CONSTITUENT S :A,B: !
    PHASE M %  1 1.0 !
    CONSTITUENT M :A,B: !
    """)

    def model(gibbs):
        class CustomModel(Model):
            contributions = [('custom', 'custom_gibbs')]

            def custom_gibbs(self, dbe):
                return gibbs([v.SiteFraction(self.phase_name, 0, s) for s in sorted(dbe.species)], v.T)
        return CustomModel
    energies = {'S': lambda x, T: (x[0] - 0.1) ** 2 + (x[1] - 0.9) ** 2 + T,
                'M': lambda x, T: (x[0] - 0.9) ** 2 + (x[1] - 0.1) ** 2 + 0.99 * T}
    models = {name: model(g)(dbf, ['A', 'B'], name) for name, g in energies.items()}
    wks = pyc.Workspace(dbf, ['A', 'B'], ['S', 'M'], {v.N: 1, v.P: 1e5, v.T: 10, v.X('A'): 0.5}, models=models)
    x_s, x_m = wks.get('X(S,A)', 'X(M,A)')
    np.testing.assert_allclose([float(x_s), float(x_m)], [0.1 - 0.03125, 0.9 - 0.03125], atol=1e-12)


def test_new_phase_search_finds_a_phase_at_the_same_composition():
    "A unary just above its melting point: the liquid (same composition, small driving force) is found."
    from pycalphad import Workspace, variables as v
    from pycalphad.core.composition_set import CompositionSet
    from pycalphad.mapping.primitives import Point
    from scipy.optimize import brentq
    from pycalphad_cg.local import global_min_candidate
    dbf = load_database('alfe.tdb')
    conds = {v.N: 1, v.P: 101325, v.T: 934}
    wks = Workspace(dbf, ['AL', 'VA'], ['FCC_A1', 'LIQUID'], conds)
    prf = wks.phase_record_factory

    def state(T):
        return np.array([float(T) if s == v.T else float(conds[s]) for s in prf.state_variables])

    def energy(name, T):
        cs = CompositionSet(prf[name])
        cs.update(np.ones(cs.phase_record.phase_dof), 1.0, state(T))
        return cs.energy
    T = brentq(lambda T: energy('LIQUID', T) - energy('FCC_A1', T), 900, 1000) + 0.05
    conds[v.T] = T
    solid = CompositionSet(prf['FCC_A1'])
    solid.update(np.ones(solid.phase_record.phase_dof), 1.0, state(T))
    point = Point(conds, np.array([solid.energy]), [], [solid])
    info = {'dbf': dbf, 'comps': ['AL', 'VA'], 'phases': ['FCC_A1', 'LIQUID'], 'models': wks.models.unwrap(),
            'phase_records': prf}
    found = global_min_candidate(point, info, tol=1e-5)
    assert found is not None and found[0].phase_record.phase_name == 'LIQUID'
    np.testing.assert_allclose(found[1], solid.energy - energy('LIQUID', T), rtol=1e-6)


def test_fallback_off_raises_in_mapping():
    "With fallback off, a problem column generation does not handle raises, although pycalphad's mapping catches it."
    import pycalphad.mapping.zpf_equilibrium as zeq
    from pycalphad import Database, Model, variables as v
    from pycalphad.mapping.starting_points import point_from_equilibrium

    class NoTemperature(Model):
        contributions = [('custom', 'custom_gibbs')]

        def custom_gibbs(self, dbe):
            return 1000 * (v.Y(self.phase_name, 0, 'AL') ** 2 + v.Y(self.phase_name, 0, 'FE') ** 2)
    dbf = Database('''
    ELEMENT AL FCC_A1 26.98 0 0 !
    ELEMENT FE BCC_A2 55.85 0 0 !
    PHASE LIQUID % 1 1 !
    CONSTITUENT LIQUID :AL,FE: !
    ''')
    conds = {v.P: 101325, v.N: 1, v.X('AL'): 0.3}
    models = {'LIQUID': NoTemperature(dbf, ['AL', 'FE'], 'LIQUID')}
    point = point_from_equilibrium(dbf, ['AL', 'FE'], ['LIQUID'], dict(conds), models=models)
    with pytest.raises(pyc.CGCapabilityError):
        with pyc.cg_mode(fallback=False):
            assert zeq.update_equilibrium_with_new_conditions(point, {**conds, v.X('AL'): 0.31}) is None
    with pytest.warns(UserWarning, match='column generation does not handle'):
        with pyc.cg_mode(fallback=True):
            assert zeq.update_equilibrium_with_new_conditions(point, {**conds, v.X('AL'): 0.31}) is not None
