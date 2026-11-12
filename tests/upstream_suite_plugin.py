"""
pytest plugin: pycalphad's own test suite with column generation for everything pycalphad_cg covers.

    cd /tmp && PYTHONPATH=/path/to/pycalphad-cg/tests python -m pytest --pyargs pycalphad.tests -p upstream_suite_plugin

Every test runs inside ``pycalphad_cg.cg_mode()``, so ``pycalphad.equilibrium``, the mapping strategies and
``binplot``/``ternplot`` use column generation, and ``pycalphad.Workspace`` is pycalphad_cg's Workspace (the test
modules import it from there). Problems column generation does not handle fall back to pycalphad's solver with a
warning; the terminal summary counts the column-generation solves and lists the fallbacks.
"""
import warnings

import pycalphad
import pytest

import pycalphad_cg
import pycalphad_cg.solve
import pycalphad_cg.workspace

_counts = {'solves': 0}
_fallbacks = []
_passes_with_cg = []

# pycalphad tests whose premise is a behaviour of pycalphad's own solver that column generation does not share. Each
# is marked xfail here; what the test is for is checked in this package's tests.
KNOWN_DIFFERENCES = {
    'test_automatic_starting_points_for_every_merged_field':
        "its premise is that a step map at 300 K steps over the narrow Au3Cu field (two fields merge into one step "
        "line); column generation resolves that field instead. The starting points and the mapped Au-rich field are "
        "checked in tests/test_mapping.py",
}


# pycalphad tests marked as expected failures of pycalphad's solver that column generation passes: run as ordinary
# tests here (an unexpected pass of a strict xfail would be reported as a failure)
PASSES_WITH_CG = {
    'test_zpf_line_from_saddle_point_configuration_reaches_invariant',   # (gh-746: L12 left on a saddle point)
}


def pytest_collection_modifyitems(config, items):
    for item in items:
        name = getattr(item, 'originalname', None) or item.name.split('[')[0]
        reason = KNOWN_DIFFERENCES.get(name)
        if reason is not None:
            item.add_marker(pytest.mark.xfail(reason=reason, strict=False))
        if name in PASSES_WITH_CG:
            item.own_markers[:] = [mark for mark in item.own_markers if mark.name != 'xfail']
            _passes_with_cg.append(item.nodeid)


class _RecordingWarnings:
    """pycalphad_cg.workspace's ``warnings``, recording the fallback warnings."""
    def __getattr__(self, name):
        return getattr(warnings, name)

    @staticmethod
    def warn(message, *args, **kwargs):
        if 'column generation does not handle' in str(message):
            _fallbacks.append(str(message))
        return warnings.warn(message, *args, **kwargs)


def pytest_configure(config):
    pycalphad.Workspace = pycalphad_cg.Workspace
    pycalphad_cg.workspace.warnings = _RecordingWarnings()
    original = pycalphad_cg.solve._solve_workspace

    def counting(*args, **kwargs):
        _counts['solves'] += 1
        return original(*args, **kwargs)
    pycalphad_cg.solve._solve_workspace = counting


@pytest.hookimpl(wrapper=True)
def pytest_runtest_call(item):
    with pycalphad_cg.cg_mode():
        return (yield)


def pytest_terminal_summary(terminalreporter):
    terminalreporter.write_line(f'pycalphad_cg {pycalphad_cg.__version__} on pycalphad {pycalphad.__version__}: '
                                f'{_counts["solves"]} column-generation solves, {len(_fallbacks)} fallbacks to '
                                f'pycalphad\'s solver')
    for reason in sorted(set(_fallbacks)):
        terminalreporter.write_line(f'    {_fallbacks.count(reason)} x {reason}')
    for nodeid in _passes_with_cg:
        terminalreporter.write_line(f'    run as an ordinary test (an expected failure of pycalphad\'s solver): {nodeid}')
