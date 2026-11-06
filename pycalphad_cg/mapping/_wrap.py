"""Column-generation versions of pycalphad's mapping strategies and compatibility functions."""
import inspect

from pycalphad_cg._cgmode import run_in_cg_mode

_CG_CLASSES = {}


def cg_strategy(cls):
    """A subclass of pycalphad strategy class ``cls`` whose methods (its own and inherited, including ``__init__``,
    except other special methods) run inside cg_mode(): its equilibria, steps and new-phase searches use column
    generation."""
    if cls in _CG_CLASSES:
        return _CG_CLASSES[cls]
    namespace = {'__doc__': f'{cls.__doc__ or ""}\n\n(pycalphad-cg: equilibria and steps by column generation)',
                 '__module__': 'pycalphad_cg.mapping'}
    seen = set()
    for klass in cls.__mro__:
        if klass is object:
            continue
        for name, attr in vars(klass).items():
            if name in seen:
                continue
            seen.add(name)
            if name.startswith('__') and name != '__init__':
                continue
            if inspect.isfunction(attr):
                namespace[name] = run_in_cg_mode(attr)
            elif isinstance(attr, (staticmethod, classmethod)) and inspect.isfunction(attr.__func__):
                namespace[name] = type(attr)(run_in_cg_mode(attr.__func__))
    new = type(cls.__name__, (cls,), namespace)
    new.__qualname__ = cls.__qualname__
    _CG_CLASSES[cls] = new
    return new


def is_strategy_class(obj):
    from pycalphad.mapping.strategy.strategy_base import MapStrategy
    return inspect.isclass(obj) and issubclass(obj, MapStrategy)


def cg_function(function):
    """``function`` (e.g. binplot) run inside cg_mode(), with the same signature and documentation."""
    return run_in_cg_mode(function)


__all__ = ['cg_strategy', 'cg_function', 'is_strategy_class']
