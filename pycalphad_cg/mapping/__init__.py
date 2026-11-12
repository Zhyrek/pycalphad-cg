"""
pycalphad's ``pycalphad.mapping``, with column generation: the strategy classes here are subclasses of pycalphad's
whose equilibria, steps and new-phase searches use column generation; everything else (plotting, primitives) is
pycalphad's own.
"""
import pycalphad.mapping as _mapping
from pycalphad.mapping.compat_api import binplot as _binplot, ternplot as _ternplot

from pycalphad_cg._aliases import public_names as _public_names, submodule as _submodule
from pycalphad_cg.mapping._wrap import cg_strategy as _cg_strategy, cg_function as _cg_function, \
    is_strategy_class as _is_strategy_class, is_deprecated as _is_deprecated, \
    deprecated_strategy as _deprecated_strategy

# (deprecated strategy classes, e.g. BinaryStrategy in pycalphad 0.12, are built on first use: see __getattr__)
globals().update({name: (_cg_strategy(value) if _is_strategy_class(value) else value)
                  for name, value in _public_names(_mapping).items()
                  if not (_is_strategy_class(value) and _is_deprecated(value))})

binplot = _cg_function(_binplot)
ternplot = _cg_function(_ternplot)


def __getattr__(name):
    cls = _deprecated_strategy(_mapping, name)
    if cls is not None:
        globals()[name] = cls
        return cls
    return _submodule(__name__, name)
