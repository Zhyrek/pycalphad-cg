"""pycalphad's ``pycalphad.mapping.strategy``, with the strategy classes in their column-generation versions (see
pycalphad_cg.mapping)."""
import pycalphad.mapping.strategy as _strategy

from pycalphad_cg._aliases import public_names as _public_names, submodule as _submodule
from pycalphad_cg.mapping._wrap import cg_strategy as _cg_strategy, is_strategy_class as _is_strategy_class

globals().update({name: (_cg_strategy(value) if _is_strategy_class(value) else value)
                  for name, value in _public_names(_strategy).items()})


def __getattr__(name):
    return _submodule(__name__, name)
