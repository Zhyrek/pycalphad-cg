"""
``pycalphad_cg.X`` for every pycalphad module ``pycalphad.X`` this package does not replace: the same module object
(not a copy), so that code written for pycalphad works with ``pycalphad`` replaced by ``pycalphad_cg`` throughout,
e.g. ``import pycalphad_cg.variables as v`` or ``from pycalphad_cg.plot.utils import phase_legend``.

The modules this package replaces are its own (``pycalphad_cg.mapping``, ``pycalphad_cg.mapping.strategy``); below
those, pycalphad's modules again (e.g. ``pycalphad_cg.mapping.plotting`` is ``pycalphad.mapping.plotting``).
"""
import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys
import types

PACKAGE = 'pycalphad_cg'
TARGET = 'pycalphad'


class _AliasLoader(importlib.abc.Loader):
    def __init__(self, target):
        self.target = target

    def create_module(self, spec):
        return importlib.import_module(self.target)

    def exec_module(self, module):
        pass   # (pycalphad's module is already executed)


class _AliasFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if not fullname.startswith(PACKAGE + '.'):
            return None
        parent_name = fullname.rpartition('.')[0]
        parent = sys.modules.get(parent_name)
        if parent is not None and getattr(parent, '__name__', '').startswith(PACKAGE):
            # a module of this package: its own submodules are loaded normally
            if importlib.machinery.PathFinder.find_spec(fullname, getattr(parent, '__path__', None)) is not None:
                return None
        target_name = TARGET + fullname[len(PACKAGE):]
        try:
            target_spec = importlib.util.find_spec(target_name)
        except (ImportError, ValueError):
            return None
        if target_spec is None:
            return None
        return importlib.machinery.ModuleSpec(fullname, _AliasLoader(target_name),
                                              is_package=target_spec.submodule_search_locations is not None)


def install():
    if not any(isinstance(finder, _AliasFinder) for finder in sys.meta_path):
        sys.meta_path.insert(0, _AliasFinder())


def public_names(module):
    """The public names of a pycalphad module that are not modules themselves (those resolve through the aliases)."""
    return {name: value for name, value in vars(module).items()
            if not name.startswith('_') and not isinstance(value, types.ModuleType)}


def submodule(package_name, name):
    """``package_name.name`` as a module (own or aliased), else AttributeError (for a module ``__getattr__``)."""
    fullname = f'{package_name}.{name}'
    if not name.startswith('__'):
        try:
            return importlib.import_module(fullname)
        except ModuleNotFoundError as err:
            if err.name != fullname:
                raise
    raise AttributeError(f'module {package_name!r} has no attribute {name!r}')


__all__ = ['install', 'public_names', 'submodule']
