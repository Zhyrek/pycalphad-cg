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
