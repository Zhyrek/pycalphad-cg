"""
pycalphad_cg's Workspace against pycalphad's own on the databases in pycalphad's test suite.

Column generation is a different minimizer, so the results are compared by energy: at every condition its molar
Gibbs energy must not be above pycalphad's (beyond rounding), and where both find the same equilibrium the chemical
potentials agree. pycalphad's grid can miss a phase, so column generation is allowed to be lower.
"""
import warnings

import numpy as np
import pytest
from pycalphad import Workspace, variables as v
from pycalphad.core.errors import ConditionError
from pycalphad.core.solver import Solver

from pycalphad_cg import CGWorkspace, CGCapabilityError, equilibrium
from pycalphad_cg.solve import _consolidate
from conftest import load_database as _db

GM_TOL = 0.01  # J/mol


def _all_phases(dbf):
    return sorted(dbf.phases.keys())


CASES = [
    ('alfe.tdb', ['AL', 'FE', 'VA'], {v.T: (500, 1901, 100), v.X('AL'): (0.02, 0.99, 0.04), v.P: 101325, v.N: 1}),
    ('cumg.tdb', ['CU', 'MG', 'VA'], {v.T: (500, 1501, 100), v.X('MG'): (0.02, 0.99, 0.04), v.P: 101325, v.N: 1}),
    ('alzn_mey.tdb', ['AL', 'ZN', 'VA'], {v.T: (400, 901, 50), v.X('ZN'): (0.02, 0.99, 0.04), v.P: 101325, v.N: 1}),
    ('FeNi_deep_branching.tdb', ['FE', 'NI', 'VA'],
     {v.T: (300, 1001, 100), v.X('NI'): (0.02, 0.99, 0.04), v.P: 101325, v.N: 1}),
    ('alcrni.tdb', ['AL', 'CR', 'NI', 'VA'],
     {v.T: 1273, v.X('AL'): (0.02, 0.5, 0.06), v.X('CR'): (0.02, 0.5, 0.06), v.P: 101325, v.N: 1}),
    ('CrFeNb_Jacob2016.tdb', ['CR', 'FE', 'NB', 'VA'],
     {v.T: 1473, v.X('CR'): (0.05, 0.6, 0.1), v.X('NB'): (0.05, 0.4, 0.1), v.P: 101325, v.N: 1}),
]


@pytest.mark.parametrize('name,comps,conds', CASES, ids=[c[0] for c in CASES])
def test_energy_not_above_pycalphad(name, comps, conds):
    dbf = _db(name)
    phases = _all_phases(dbf)
    ref = Workspace(dbf, comps, phases, conds).eq
    cg = CGWorkspace(dbf, comps, phases, conds).eq
    gm_ref, gm_cg = np.squeeze(ref.GM), np.squeeze(cg.GM)
    valid = ~np.isnan(gm_ref)
    assert np.all(~np.isnan(gm_cg[valid]))
    assert np.all(gm_cg[valid] <= gm_ref[valid] + GM_TOL), np.max(gm_cg[valid] - gm_ref[valid])
    # where the two agree on the energy, they agree on the chemical potentials
    same = valid & (np.abs(gm_cg - gm_ref) < GM_TOL)
    mu_ref, mu_cg = np.squeeze(ref.MU)[same], np.squeeze(cg.MU)[same]
    np.testing.assert_allclose(mu_cg, mu_ref, atol=1.0)


def test_equilibrium_wrapper_and_properties():
    dbf = _db('alfe.tdb')
    phases = _all_phases(dbf)
    conds = {v.T: (800, 1201, 200), v.X('AL'): (0.1, 0.6, 0.1), v.P: 101325, v.N: 1}
    from pycalphad import equilibrium as pycalphad_equilibrium
    ref = pycalphad_equilibrium(dbf, ['AL', 'FE', 'VA'], phases, conds, output='HM')
    cg = equilibrium(dbf, ['AL', 'FE', 'VA'], phases, conds, output='HM')
    assert set(ref.data_vars) == set(cg.data_vars)
    assert ref.GM.dims == cg.GM.dims
    np.testing.assert_allclose(cg.HM.values, ref.HM.values, atol=1.0)


def test_workspace_get():
    dbf = _db('cumg.tdb')
    conds = {v.T: 1000, v.X('MG'): 0.3, v.P: 101325, v.N: 1}
    phases = _all_phases(dbf)
    ref = Workspace(dbf, ['CU', 'MG', 'VA'], phases, conds)
    cg = CGWorkspace(dbf, ['CU', 'MG', 'VA'], phases, conds)
    for prop in ['GM', 'HM', 'MU(MG)', 'NP(*)', 'X(LIQUID,MG)']:
        np.testing.assert_allclose(np.nan_to_num(cg.get(prop)), np.nan_to_num(ref.get(prop)), atol=1e-3, rtol=1e-6)


def test_input_errors_match_pycalphad():
    dbf = _db('alfe.tdb')
    phases = _all_phases(dbf)
    # missing state variable
    with pytest.raises(ConditionError):
        CGWorkspace(dbf, ['AL', 'FE', 'VA'], phases, {v.X('AL'): 0.3, v.P: 101325, v.N: 1}).eq
    # N != 1
    with pytest.raises(ConditionError):
        CGWorkspace(dbf, ['AL', 'FE', 'VA'], phases, {v.T: 1000, v.X('AL'): 0.3, v.P: 101325, v.N: 2}).eq
    # wrong number of composition conditions
    with pytest.raises(ValueError, match='degrees of freedom'):
        CGWorkspace(dbf, ['AL', 'FE', 'VA'], phases, {v.T: 1000, v.P: 101325, v.N: 1}).eq


class _CustomSolver(Solver):
    pass


def test_fallback_to_pycalphad():
    dbf = _db('alfe.tdb')
    phases = _all_phases(dbf)
    conds = {v.T: 1000, v.X('AL'): 0.3, v.P: 101325, v.N: 1}
    with pytest.warns(UserWarning, match='column generation does not handle'):
        gm = CGWorkspace(dbf, ['AL', 'FE', 'VA'], phases, conds, solver=_CustomSolver()).get('GM')
    np.testing.assert_allclose(gm, Workspace(dbf, ['AL', 'FE', 'VA'], phases, conds).get('GM'))
    with pytest.raises(CGCapabilityError):
        CGWorkspace(dbf, ['AL', 'FE', 'VA'], phases, conds, solver=_CustomSolver(), cg_fallback=False).eq


def test_threads_give_the_same_result(monkeypatch):
    "Several worker threads (verified) give the single-thread result."
    import pycalphad_cg.solve as solve
    created = []

    class CountingWorker(solve._Worker):
        def __init__(self, index, *args, **kwargs):
            created.append(index)
            super().__init__(index, *args, **kwargs)
    monkeypatch.setattr(solve, '_Worker', CountingWorker)
    dbf = _db('alcrni.tdb')
    phases = _all_phases(dbf)
    conds = {v.T: 1273, v.X('AL'): (0.05, 0.5, 0.05), v.X('CR'): (0.05, 0.5, 0.05), v.P: 101325, v.N: 1}
    one = CGWorkspace(dbf, ['AL', 'CR', 'NI', 'VA'], phases, conds, cg_threads=1).get('GM')
    assert created == [0]
    created.clear()
    four = CGWorkspace(dbf, ['AL', 'CR', 'NI', 'VA'], phases, conds, cg_threads=4).get('GM')
    assert created == [0, 1, 2, 3]
    valid = ~np.isnan(np.asarray(Workspace(dbf, ['AL', 'CR', 'NI', 'VA'], phases, conds).get('GM')))
    assert valid.sum() > 50 and np.all(np.isfinite(np.asarray(one)[valid]))
    np.testing.assert_allclose(four, one, atol=1e-6)


def test_fixed_chemical_potential_and_mass_fraction():
    dbf = _db('alfe.tdb')
    phases = _all_phases(dbf)
    comps = ['AL', 'FE', 'VA']
    base = {v.T: 1200, v.P: 101325, v.N: 1}
    for conds in ({**base, v.MU('AL'): -90000}, {**base, v.W('AL'): 0.2}):
        ref = Workspace(dbf, comps, phases, conds).get('GM', 'X(AL)')
        cg = CGWorkspace(dbf, comps, phases, conds).get('GM', 'X(AL)')
        np.testing.assert_allclose(cg, ref, atol=1e-3, rtol=1e-6)


def test_consolidate_keeps_mass_balance_and_energy():
    "Merging two points of one phase within COMP_DIFFERENCE_TOL keeps the overall composition and GM."
    mu = np.array([-30000.0, -60000.0])
    XA, XB, XC = np.array([0.90, 0.10]), np.array([0.90004, 0.09996]), np.array([0.3, 0.7])
    result = ([0, 0, 1], [np.array([0.9, 0.1]), np.array([0.90004, 0.09996]), np.array([0.3, 0.7])],
              [XA, XB, XC], [0.4, 0.3, 0.3], [mu @ XA, mu @ XB, mu @ XC])
    z = 0.4 * XA + 0.3 * XB + 0.3 * XC
    gm = 0.4 * mu @ XA + 0.3 * mu @ XB + 0.3 * mu @ XC
    phs, Y, X, n, G = _consolidate(result)
    assert phs == [0, 1]
    np.testing.assert_allclose(sum(a * x for a, x in zip(n, X)), z, atol=1e-15)
    assert abs(sum(a * g for a, g in zip(n, G)) - gm) < 1e-9


def _q_database(n):
    "A liquid with an associate species Q = AL1NI{n}: the same energy expressions, different formula moles."
    from pycalphad import Database
    return Database(f'''
    ELEMENT AL FCC_A1 26.98 0 0 !
    ELEMENT NI FCC_A1 58.69 0 0 !
    SPECIES Q AL1NI{n} !
    TYPE_DEFINITION % SEQ * !
    PHASE LIQUID % 1 1 !
    CONSTITUENT LIQUID :AL,NI,Q: !
    PARAMETER G(LIQUID,AL;0) 1 0; 6000 N !
    PARAMETER G(LIQUID,NI;0) 1 0; 6000 N !
    PARAMETER G(LIQUID,Q;0) 1 -10000; 6000 N !
    ''')


def test_phase_cache_tells_apart_species_with_other_constituents():
    "Compiled phases are reused across calls only for the same structure (here: what the species Q is made of)."
    conds = {v.N: 1, v.P: 101325, v.T: 300, v.X('AL'): 0.4}
    for n in (1, 2):
        dbf = _q_database(n)
        cg = CGWorkspace(dbf, ['AL', 'NI'], ['LIQUID'], conds).get('GM')
        ref = Workspace(dbf, ['AL', 'NI'], ['LIQUID'], conds).get('GM')
        np.testing.assert_allclose(cg, ref, atol=1e-3)


@pytest.mark.parametrize('condition', [(v.X('AL'), 1 - 1e-10), (v.X('FE'), 1e-10), (v.W('AL'), 1e-10),
                                       (v.W('AL'), 1 - 1e-10)], ids=['X_AL', 'X_FE', 'W_AL', 'W_AL_high'])
def test_trace_compositions_in_any_form(condition):
    "Trace components implied by any condition (the dependent component, mass fractions) are solved accurately."
    dbf = _db('alfe.tdb')
    conds = {v.N: 1, v.T: 1200, v.P: 101325, condition[0]: condition[1]}
    ref = Workspace(dbf, ['AL', 'FE', 'VA'], ['LIQUID'], conds).eq
    cg = CGWorkspace(dbf, ['AL', 'FE', 'VA'], ['LIQUID'], conds).eq
    np.testing.assert_allclose(cg.X.ravel()[:2], ref.X.ravel()[:2], rtol=1e-4)
    np.testing.assert_allclose(cg.MU.ravel(), ref.MU.ravel(), atol=0.1)


@pytest.mark.parametrize('y_mg', [0.1, 0.3])
def test_phase_local_conditions_give_pycalphads_potentials(y_mg):
    "Where phase-local conditions leave the chemical potentials undetermined, both give the minimum-norm ones."
    dbf = _db('cumg.tdb')
    conds = {v.N: 1, v.P: 101325, v.T: 1080, v.X('MG'): y_mg, v.Y('LIQUID', 0, 'MG'): y_mg}
    ref = Workspace(dbf, ['CU', 'MG'], ['LIQUID'], conds).eq
    cg = CGWorkspace(dbf, ['CU', 'MG'], ['LIQUID'], conds).eq
    np.testing.assert_allclose(cg.MU.ravel(), ref.MU.ravel(), atol=1e-3)


def test_system_amount_in_other_units():
    dbf = _db('alfe.tdb')
    conds = {v.N['mmol']: 1000, v.T: 1200, v.P: 101325, v.X('AL'): 0.2}
    ref = Workspace(dbf, ['AL', 'FE', 'VA'], ['LIQUID'], conds).get('GM')
    np.testing.assert_allclose(CGWorkspace(dbf, ['AL', 'FE', 'VA'], ['LIQUID'], conds).get('GM'), ref, atol=1e-3)


def test_configured_solver_is_honoured():
    "A pycalphad Solver with non-default settings is the caller's choice: pycalphad's solver is used with it."
    dbf = _db('alfe.tdb')
    conds = {v.T: 1000, v.X('AL'): 0.3, v.P: 101325, v.N: 1}
    with pytest.warns(UserWarning, match='column generation does not handle'):
        CGWorkspace(dbf, ['AL', 'FE', 'VA'], _all_phases(dbf), conds, solver=Solver(remove_metastable=False)).eq
    with warnings.catch_warnings():
        warnings.simplefilter('error')   # (a default-configured Solver is replaced without a warning)
        CGWorkspace(dbf, ['AL', 'FE', 'VA'], _all_phases(dbf), conds, solver=Solver(verbose=False)).eq


def test_energies_far_from_zero():
    "The linear program does not depend on the energy reference (e.g. all end members shifted by 1e7 J/mol)."
    from pycalphad import Database
    template = '''
    ELEMENT AL FCC_A1 26.98 0 0 !
    ELEMENT NI FCC_A1 58.69 0 0 !
    TYPE_DEFINITION % SEQ * !
    PHASE LIQUID % 1 1 !
    CONSTITUENT LIQUID :AL,NI: !
    PARAMETER G(LIQUID,AL;0) 1 {0}; 6000 N !
    PARAMETER G(LIQUID,NI;0) 1 {0}; 6000 N !
    '''
    conds = {v.N: 1, v.P: 101325, v.T: 300, v.X('AL'): 0.5}
    gm = [float(CGWorkspace(Database(template.format(shift)), ['AL', 'NI'], ['LIQUID'], conds).get('GM'))
          for shift in (0, 10000000)]
    np.testing.assert_allclose(gm[1] - gm[0], 1e7, rtol=1e-12)


def test_ratio_condition_is_met():
    "A ratio of mole fractions as a condition (the denominator taken by element, not by position)."
    from pycalphad import Database
    dbf = Database('''
    ELEMENT AL FCC_A1 26.98 0 0 !
    ELEMENT NI FCC_A1 58.69 0 0 !
    ELEMENT ZN HCP_A3 65.38 0 0 !
    TYPE_DEFINITION % SEQ * !
    PHASE LIQUID % 1 1 !
    CONSTITUENT LIQUID :AL,NI,ZN: !
    ''')
    conds = {v.N: 1, v.P: 101325, v.T: 300, v.X('AL'): 0.2, v.X('NI') / v.X('ZN'): 2}
    x = CGWorkspace(dbf, ['AL', 'NI', 'ZN'], ['LIQUID'], conds).eq.X.reshape(-1, 3)[0]
    np.testing.assert_allclose(x, [0.2, 0.8 * 2 / 3, 0.8 / 3], atol=1e-9)
