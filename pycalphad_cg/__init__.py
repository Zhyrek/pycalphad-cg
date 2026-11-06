"""
pycalphad with column generation.

    import pycalphad_cg as pyc       # in place of: import pycalphad as pyc
    from pycalphad_cg import Database, equilibrium, binplot, variables as v

This module offers pycalphad's API. Equilibria come from column generation: ``equilibrium``, ``Workspace``, the
mapping strategies (``pycalphad_cg.mapping``), ``binplot`` and ``ternplot``. Everything else is pycalphad's own
object (``pycalphad_cg.Database is pycalphad.Database``), and so are the modules: ``pycalphad_cg.X`` is
``pycalphad.X`` for every module this package does not replace.

Additional names: ``options`` (package-wide defaults), ``cg_mode`` (column generation for pycalphad's own functions
within a ``with`` block), ``CGWorkspace`` (= ``Workspace``), ``CGCapabilityError``.
"""
from pycalphad_cg._version import __version__
from pycalphad_cg import _aliases

_aliases.install()

import pycalphad as _pycalphad  # noqa: E402  (pycalphad imports this package as a plugin; see _aliases)
from pycalphad import (Database, Model, ReferenceState, calculate, eqplot, as_property,  # noqa: E402
                       variables, CalculateError, ConditionError, DofError, EquilibriumError)
from pycalphad_cg._options import options  # noqa: E402
from pycalphad_cg._cgmode import cg_mode  # noqa: E402
from pycalphad_cg.solve import CGCapabilityError  # noqa: E402
from pycalphad_cg.workspace import Workspace, CGWorkspace, equilibrium  # noqa: E402
from pycalphad_cg import mapping  # noqa: E402
from pycalphad_cg.mapping import binplot, ternplot  # noqa: E402

v = variables
pycalphad_version = _pycalphad.__version__

__all__ = ['Database', 'Model', 'ReferenceState', 'calculate', 'equilibrium', 'Workspace', 'eqplot', 'as_property',
           'binplot', 'ternplot', 'variables', 'v', 'CalculateError', 'ConditionError', 'DofError',
           'EquilibriumError', 'mapping', 'options', 'cg_mode', 'CGWorkspace', 'CGCapabilityError', '__version__',
           'pycalphad_version']


def __getattr__(name):
    # names pycalphad gains in later versions, and its modules (see _aliases)
    if not name.startswith('_'):
        value = getattr(_pycalphad, name, None)
        if value is not None and not isinstance(value, type(_pycalphad)):
            return value
    return _aliases.submodule(__name__, name)
