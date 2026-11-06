"""pycalphad's ``pycalphad.mapping.strategy``, with the strategy classes in their column-generation versions (see
pycalphad_cg.mapping)."""
import pycalphad.mapping.strategy as _strategy

from pycalphad_cg._aliases import public_names as _public_names, submodule as _submodule
from pycalphad_cg.mapping._wrap import all_cg_strategies as _all_cg_strategies

globals().update(_public_names(_strategy))
globals().update(_all_cg_strategies())   # (every column-generation strategy class lives here, see cg_strategy)


def __getattr__(name):
    return _submodule(__name__, name)
