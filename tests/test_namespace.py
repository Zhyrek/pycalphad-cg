"""pycalphad_cg as a drop-in for the pycalphad namespace."""
import subprocess
import sys
import tempfile
import threading

import numpy as np
import pycalphad
import pycalphad.core.equilibrium
import pycalphad.mapping.starting_points
import pycalphad.mapping.zpf_equilibrium
import pytest

import pycalphad_cg as pyc
from conftest import load_database


def test_pycalphad_objects_are_shared():
    for name in ['Database', 'Model', 'ReferenceState', 'calculate', 'eqplot', 'as_property', 'variables',
                 'ConditionError', 'CalculateError', 'DofError', 'EquilibriumError']:
        assert getattr(pyc, name) is getattr(pycalphad, name), name
    assert pyc.v is pycalphad.variables


def test_column_generation_entry_points():
    assert issubclass(pyc.Workspace, pycalphad.Workspace) and pyc.Workspace is not pycalphad.Workspace
    assert pyc.CGWorkspace is pyc.Workspace
    assert pyc.equilibrium is not pycalphad.equilibrium
    assert pyc.binplot is not pycalphad.binplot and pyc.ternplot is not pycalphad.ternplot
    for name in ['BinaryStrategy', 'TernaryStrategy', 'IsoplethStrategy', 'StepStrategy']:
        if hasattr(pycalphad.mapping, name):
            cg_class = getattr(pyc.mapping, name)
            assert cg_class.__mro__[1] is getattr(pycalphad.mapping, name)
            assert getattr(pyc.mapping.strategy, name) is cg_class


def test_modules_alias_pycalphad():
    import pycalphad_cg.variables as variables
    import pycalphad_cg.core.workspace
    from pycalphad_cg.plot.utils import phase_legend
    from pycalphad.plot.utils import phase_legend as pycalphad_phase_legend
    import pycalphad.mapping.plotting
    assert variables is pycalphad.variables
    assert pycalphad_cg.core.workspace is sys.modules['pycalphad.core.workspace']
    assert phase_legend is pycalphad_phase_legend
    assert pyc.mapping.plotting is pycalphad.mapping.plotting
    assert pyc.plot is pycalphad.plot
    assert pyc.mapping.plot_binary is pycalphad.mapping.plot_binary


def test_missing_names_raise():
    with pytest.raises(AttributeError):
        pyc.no_such_name
    with pytest.raises(ModuleNotFoundError):
        import pycalphad_cg.no_such_module  # noqa: F401
    assert not hasattr(pyc, '__no_such_dunder__')


@pytest.mark.parametrize('first', ['pycalphad', 'pycalphad_cg'])
def test_import_order(first):
    "Either import order works (pycalphad imports pycalphad_* packages as plugins), and plain pycalphad is unchanged."
    code = (f'import {first}; import pycalphad, pycalphad_cg, pycalphad.core.equilibrium as e; '
            'assert e.Workspace is pycalphad.Workspace and pycalphad.equilibrium is e.equilibrium; '
            'assert pycalphad_cg.equilibrium is not pycalphad.equilibrium; print("ok")')
    # (in a neutral directory: python -c puts the working directory first on the path)
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=tempfile.gettempdir())
    assert out.returncode == 0 and out.stdout.strip() == 'ok', out.stderr


def _hooked_names():
    return [pycalphad.core.equilibrium.Workspace, pycalphad.mapping.starting_points.Workspace,
            pycalphad.mapping.zpf_equilibrium.Solver, pycalphad.mapping.zpf_equilibrium._find_global_min_cs]


def test_cg_mode_is_per_thread_and_restores(cg_calls):
    dbf = load_database('alfe.tdb')
    args = (dbf, ['AL', 'FE', 'VA'], list(dbf.phases), {pyc.v.T: 1000, pyc.v.X('AL'): 0.3, pyc.v.P: 101325,
                                                         pyc.v.N: 1})
    before = _hooked_names()
    inside, release = threading.Event(), threading.Event()
    results = {}

    def cg_thread():
        with pyc.cg_mode():
            results['cg'] = pycalphad.equilibrium(*args)   # pycalphad's own function, in cg mode
            inside.set()
            release.wait(60)

    worker = threading.Thread(target=cg_thread)
    worker.start()
    assert inside.wait(120)
    assert _hooked_names() != before                       # the names are redirected while a thread is in cg mode
    assert cg_calls['workspace'] == 1
    results['plain'] = pycalphad.equilibrium(*args)        # another thread meanwhile: pycalphad's own solver
    assert cg_calls['workspace'] == 1
    release.set()
    worker.join()
    assert _hooked_names() == before                       # and restored when no thread is
    np.testing.assert_allclose(results['cg'].GM.values, results['plain'].GM.values, atol=1e-3)


def test_settings_nest(monkeypatch):
    from pycalphad_cg import _cgmode
    with pyc.cg_mode(threads=2):
        assert _cgmode._settings() == (2, pyc.options.fallback)
        with pyc.cg_mode(fallback=False):
            assert _cgmode._settings() == (2, False)
        assert _cgmode._settings() == (2, pyc.options.fallback)
    monkeypatch.setattr(pyc.options, 'threads', 3)
    with pyc.cg_mode():
        assert _cgmode._settings()[0] == 3
    assert not _cgmode.in_cg_mode()
