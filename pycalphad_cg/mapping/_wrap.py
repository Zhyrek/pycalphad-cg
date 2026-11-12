"""Column-generation versions of pycalphad's mapping strategies and compatibility functions."""
import inspect
import warnings

from pycalphad_cg._cgmode import run_in_cg_mode

_CG_CLASSES = {}


def cg_strategy(cls):
    """A subclass of pycalphad strategy class ``cls`` whose methods (its own and inherited, including ``__init__``,
    except other special methods) run inside cg_mode(): its equilibria, steps and new-phase searches use column
    generation. The column-generation versions keep pycalphad's hierarchy (the version of a subclass is a subclass
    of the version of its base), and they live in pycalphad_cg.mapping.strategy (for pickling)."""
    if cls in _CG_CLASSES:
        return _CG_CLASSES[cls]
    namespace = {'__doc__': f'{cls.__doc__ or ""}\n\n(pycalphad-cg: equilibria and steps by column generation)',
                 '__module__': 'pycalphad_cg.mapping.strategy'}
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
    bases = (cls,) + tuple(cg_strategy(base) for base in cls.__bases__ if is_strategy_class(base))
    new = type(cls.__name__, bases, namespace)
    new.__qualname__ = cls.__name__
    _CG_CLASSES[cls] = new
    return new


def is_deprecated(cls):
    """Whether pycalphad marks the class deprecated (``warnings.deprecated``, e.g. BinaryStrategy in pycalphad 0.12),
    which warns when the class is subclassed, not only when it is used."""
    return getattr(cls, '__deprecated__', None) is not None


def cg_strategy_quietly(cls):
    """cg_strategy without the deprecation warning that subclassing a deprecated class gives: the column-generation
    version is built on request, and using it still warns, as pycalphad's class does."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        return cg_strategy(cls)


def all_cg_strategies(include_deprecated=False):
    """The column-generation versions of pycalphad's strategy classes (by name); deprecated ones only on request
    (importing this package must not warn: pycalphad imports it as a plugin)."""
    import pycalphad.mapping as mapping
    import pycalphad.mapping.strategy as strategy
    out = {}
    for module in (mapping, strategy):
        for name, value in vars(module).items():
            if not name.startswith('_') and is_strategy_class(value) and (include_deprecated or not is_deprecated(value)):
                out[name] = cg_strategy(value)
    return out


def deprecated_strategy(module, name):
    """The column-generation version of the deprecated strategy class ``module.name`` of pycalphad, or None."""
    value = getattr(module, name, None)
    if is_strategy_class(value) and is_deprecated(value):
        return cg_strategy_quietly(value)
    return None


def is_strategy_class(obj):
    from pycalphad.mapping.strategy.strategy_base import MapStrategy
    return inspect.isclass(obj) and issubclass(obj, MapStrategy)


def cg_function(function):
    """``function`` (e.g. binplot) run inside cg_mode(), with the same signature and documentation."""
    return run_in_cg_mode(function)


__all__ = ['cg_strategy', 'cg_function', 'is_strategy_class', 'all_cg_strategies', 'deprecated_strategy']
