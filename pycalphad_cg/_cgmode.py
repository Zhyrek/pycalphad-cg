"""
Column generation inside pycalphad's own code paths.

pycalphad has no hook for an outside solver in ``equilibrium`` or in mapping: ``pycalphad.equilibrium`` builds a
``Workspace``, the mapping strategies' starting points build a ``Workspace``, and their stepping builds pycalphad's
``Solver`` and searches a grid for new phases (``_find_global_min_cs``), each by a module-level name. While a thread
is inside ``cg_mode()``, those four names resolve to the column-generation versions for that thread only: the names
are replaced by dispatchers while any thread is in the mode (and restored afterwards), and the dispatchers hand
other threads pycalphad's own objects. Nothing changes outside the mode, also not when pycalphad imports this package
as a plugin.
"""
import contextlib
import functools
import inspect
import threading
import warnings

from pycalphad_cg._options import options

__all__ = ['cg_mode', 'in_cg_mode']

_local = threading.local()
_lock = threading.Lock()
_users = 0          # threads (and nesting levels) inside cg_mode
_installed = {}     # (module, name) -> pycalphad's original object
_warned = set()


def in_cg_mode():
    """Whether the calling thread is inside cg_mode()."""
    return getattr(_local, 'settings', None) is not None


def _settings():
    threads, fallback = _local.settings
    return (options.threads if threads is None else threads), (options.fallback if fallback is None else fallback)


def _cg_workspace(original, *args, **kwargs):
    from pycalphad_cg.workspace import Workspace
    threads, fallback = _settings()
    return Workspace(*args, cg_threads=threads, cg_fallback=fallback, **kwargs)


def _cg_solver(original, *args, **kwargs):
    from pycalphad_cg.local import CGSolver
    return CGSolver(*args, **kwargs)


def _cg_global_min(original, *args, **kwargs):
    from pycalphad_cg.local import global_min_candidate
    from pycalphad_cg.solve import CGCapabilityError
    bound = inspect.signature(original).bind(*args, **kwargs)
    bound.apply_defaults()
    arguments = bound.arguments
    try:
        return global_min_candidate(arguments['point'], arguments['system_info'], arguments['tol'])
    except CGCapabilityError:
        return original(*args, **kwargs)


def _targets():
    import pycalphad.core.equilibrium as equilibrium_module
    import pycalphad.mapping.starting_points as starting_points
    import pycalphad.mapping.zpf_equilibrium as zpf_equilibrium
    return [(equilibrium_module, 'Workspace', _cg_workspace),
            (starting_points, 'Workspace', _cg_workspace),
            (zpf_equilibrium, 'Solver', _cg_solver),
            (zpf_equilibrium, '_find_global_min_cs', _cg_global_min)]


def _dispatcher(original, cg_version):
    @functools.wraps(original)
    def dispatch(*args, **kwargs):
        if in_cg_mode():
            return cg_version(original, *args, **kwargs)
        return original(*args, **kwargs)
    dispatch._pycalphad_cg_original = original
    return dispatch


def _install():
    for module, name, cg_version in _targets():
        current = getattr(module, name, None)
        if current is None:
            key = f'{module.__name__}.{name}'
            if key not in _warned:
                _warned.add(key)
                warnings.warn(f'pycalphad-cg: {key} not found in this pycalphad version; that path keeps using '
                              f'pycalphad\'s own solver', RuntimeWarning, stacklevel=4)
            continue
        if hasattr(current, '_pycalphad_cg_original'):
            continue
        _installed[(module, name)] = current
        setattr(module, name, _dispatcher(current, cg_version))


def _uninstall():
    for (module, name), original in _installed.items():
        if getattr(getattr(module, name, None), '_pycalphad_cg_original', None) is original:
            setattr(module, name, original)
    _installed.clear()


@contextlib.contextmanager
def cg_mode(threads=None, fallback=None):
    """
    Inside this context, pycalphad's ``equilibrium`` and mapping (strategies, ``binplot``, ``ternplot``) called from
    this thread use column generation. ``threads`` and ``fallback`` default to ``pycalphad_cg.options`` (an outer
    cg_mode's values, when nested).
    """
    global _users
    previous = getattr(_local, 'settings', None)
    if previous is not None:
        threads = previous[0] if threads is None else threads
        fallback = previous[1] if fallback is None else fallback
    with _lock:
        if _users == 0:
            _install()
        _users += 1
    _local.settings = (threads, fallback)
    try:
        yield
    finally:
        _local.settings = previous
        with _lock:
            _users -= 1
            if _users == 0:
                _uninstall()


def run_in_cg_mode(function):
    """``function`` with every call inside cg_mode()."""
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        with cg_mode():
            return function(*args, **kwargs)
    return wrapper
