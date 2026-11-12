"""pycalphad's ``pycalphad.mapping.strategy``, with the strategy classes in their column-generation versions (see
pycalphad_cg.mapping)."""
import pycalphad.mapping.strategy as _strategy

from pycalphad_cg._aliases import public_names as _public_names, submodule as _submodule
from pycalphad_cg.mapping._wrap import all_cg_strategies as _all_cg_strategies, \
    deprecated_strategy as _deprecated_strategy, is_deprecated as _is_deprecated, is_strategy_class as _is_strategy_class

globals().update({name: value for name, value in _public_names(_strategy).items()
                  if not (_is_strategy_class(value) and _is_deprecated(value))})
# (every column-generation strategy class lives here, see cg_strategy; deprecated ones are built on first use)
globals().update(_all_cg_strategies())


def __getattr__(name):
    import pycalphad.mapping
    for module in (_strategy, pycalphad.mapping):
        cls = _deprecated_strategy(module, name)
        if cls is not None:
            globals()[name] = cls
            return cls
    return _submodule(__name__, name)
