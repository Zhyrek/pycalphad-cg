from importlib.resources import files

import pytest
from pycalphad import Database

DATABASES = files('pycalphad.tests.databases')


def load_database(name):
    return Database(str(DATABASES / name))


@pytest.fixture
def cg_calls(monkeypatch):
    """Counts of the column-generation entry points reached: Workspace solves, local (mapping) solves and new-phase
    searches."""
    import pycalphad_cg.local as local
    import pycalphad_cg.solve as solve
    calls = {'workspace': 0, 'local': 0, 'global_min': 0}

    def counting(key, function):
        def wrapper(*args, **kwargs):
            calls[key] += 1
            return function(*args, **kwargs)
        return wrapper
    monkeypatch.setattr(solve, '_solve_workspace', counting('workspace', solve._solve_workspace))
    monkeypatch.setattr(local.CGSolver, '_solve', counting('local', local.CGSolver._solve))
    monkeypatch.setattr(local, '_global_min_candidate', counting('global_min', local._global_min_candidate))
    return calls
