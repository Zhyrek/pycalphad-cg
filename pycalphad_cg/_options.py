"""Package-wide defaults (``pycalphad_cg.options``)."""


class Options:
    """
    Defaults for every column-generation calculation; a Workspace's ``cg_threads`` / ``cg_fallback`` (or the same
    keyword arguments of ``equilibrium``) override them for that calculation.

    threads : int or None
        Worker threads per calculation. None: the ``PYCALPHAD_CG_THREADS`` environment variable, else all cores.
    fallback : bool
        True: problems column generation does not handle are solved by pycalphad's own solver, with a warning.
        False: they raise ``CGCapabilityError``.
    """
    def __init__(self):
        self.threads = None
        self.fallback = True

    def __repr__(self):
        return f'Options(threads={self.threads!r}, fallback={self.fallback!r})'


options = Options()
