"""
pycalphad's ``pycalphad.mapping``, with column generation: the strategy classes here are subclasses of pycalphad's
whose equilibria, steps and new-phase searches use column generation; everything else (plotting, primitives) is
pycalphad's own.
"""
import pycalphad.mapping as _mapping
from pycalphad.mapping.compat_api import binplot as _binplot, ternplot as _ternplot

from pycalphad_cg._aliases import public_names as _public_names, submodule as _submodule
from pycalphad_cg.mapping._wrap import cg_strategy as _cg_strategy, cg_function as _cg_function, \
    is_strategy_class as _is_strategy_class

globals().update({name: (_cg_strategy(value) if _is_strategy_class(value) else value)
                  for name, value in _public_names(_mapping).items()})

binplot = _cg_function(_binplot)
ternplot = _cg_function(_ternplot)


def __getattr__(name):
    return _submodule(__name__, name)
