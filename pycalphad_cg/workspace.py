"""
pycalphad's Workspace with its equilibria from column generation, and ``equilibrium`` on it.

pycalphad's own Workspace finds equilibria by sampling every phase on a grid, taking the lower convex hull of the
samples as starting points and refining them (``Workspace.recompute``). This Workspace overrides ``recompute`` only:
it returns a dataset in the same layout, so everything that reads ``wks.eq`` (``wks.get(...)``, properties, Jansson
derivatives, reference states) is pycalphad's, unchanged.
"""
import functools
import warnings

import numpy as np
from pycalphad.codegen.phase_record_factory import PhaseRecordFactory
from pycalphad.core.equilibrium import equilibrium as _pycalphad_equilibrium
from pycalphad.core.solver import Solver
from pycalphad.core.workspace import Workspace as _PycalphadWorkspace

from pycalphad_cg._cgmode import cg_mode
from pycalphad_cg._options import options
from pycalphad_cg.solve import solve_workspace, CGCapabilityError

__all__ = ['Workspace', 'CGWorkspace', 'equilibrium']


class Workspace(_PycalphadWorkspace):
    """
    pycalphad's Workspace, with the equilibria computed by column generation (see the pycalphad-cg README).

    Takes the arguments of ``pycalphad.Workspace``, and two more:

    cg_threads : int or None
        Worker threads for the conditions (default: ``pycalphad_cg.options.threads``).
    cg_fallback : bool or None
        Whether problems column generation does not handle (see the README) are solved by pycalphad's own
        ``Workspace.recompute``, with a warning, or raise CGCapabilityError (default:
        ``pycalphad_cg.options.fallback``).
    """
    def __init__(self, *args, cg_threads=None, cg_fallback=None, **kwargs):
        # (plain attributes: changing them does not change the equilibrium, so it is not recomputed)
        self.cg_threads = cg_threads
        self.cg_fallback = cg_fallback
        super().__init__(*args, **kwargs)

    def recompute(self):
        self._validate_inputs()
        _check_degrees_of_freedom(self.conditions, self.phase_record_factory)
        fallback = options.fallback if self.cg_fallback is None else self.cg_fallback
        try:
            self._check_supported()
            self.phase_record_factory.update_parameters(self.parameters.unwrap())
            threads = options.threads if self.cg_threads is None else self.cg_threads
            return solve_workspace(self, threads=threads)
        except CGCapabilityError as err:
            if not fallback:
                raise
            warnings.warn(f'column generation does not handle this problem ({err}); using pycalphad\'s solver',
                          stacklevel=2)
            return super().recompute()

    def _validate_inputs(self):
        """The input checks pycalphad makes in ``calculate`` (which column generation does not use), with the same
        errors."""
        from pycalphad.core.errors import ConditionError
        from pycalphad.core.utils import filter_phases, unpack_species
        comps = sorted(unpack_species(self.database, self.components))
        possible = filter_phases(self.database, comps)
        if len(possible) == 0:
            raise ConditionError(f'There are no phases in the Database that can be active with components {comps}')
        if len(filter_phases(self.database, comps, self.phases)) == 0:
            raise ConditionError(f'None of the passed phases ({self.phases}) are active. List of possible phases: '
                                 f'{possible}.')
        conds = {str(key): value for key, value in self.conditions.items()}
        if 'N' in conds and np.any(np.asarray(conds['N']) != 1):
            raise ConditionError(f'N!=1 is not yet supported, got N={conds["N"]}')
        missing = sorted({str(sv) for sv in self.phase_record_factory.state_variables} - set(conds) - {'N'})
        if missing:
            raise ConditionError(f'The following state variable(s) are required by the Model(s) for the active '
                                 f'phases but were not specified: {missing}.')

    def _check_supported(self):
        if self.calc_opts and set(self.calc_opts) - {'pdens'}:
            raise CGCapabilityError(f'calc_opts {sorted(set(self.calc_opts) - {"pdens"})}')
        if type(self.solver) is not Solver:
            raise CGCapabilityError('a custom solver object (pycalphad refines its own starting points with it)')
        if type(self.phase_record_factory) is not PhaseRecordFactory:
            raise CGCapabilityError('custom phase records')


CGWorkspace = Workspace


def _check_degrees_of_freedom(conditions, phase_record_factory):
    """The check pycalphad makes in ``starting_point`` (which column generation does not use): the conditions on
    components must leave exactly one dependent component."""
    from pycalphad.property_framework.computed_property import LinearCombination
    elements = phase_record_factory.nonvacant_elements
    number_dof = len(elements) - 1
    for key in conditions.keys():
        if not (hasattr(key, 'species') or isinstance(key, LinearCombination)):
            continue
        if hasattr(key, 'species') and getattr(key, 'phase_name', None) is not None:
            continue   # phase-local conditions do not reduce the degrees of freedom
        number_dof -= 1
    if number_dof != 0:
        raise ValueError('Number of degrees of freedom is not zero')


@functools.wraps(_pycalphad_equilibrium)
def equilibrium(*args, cg_threads=None, cg_fallback=None, **kwargs):
    with cg_mode(cg_threads, cg_fallback):
        return _pycalphad_equilibrium(*args, **kwargs)


equilibrium.__doc__ = ("pycalphad.equilibrium with the equilibria computed by column generation; same arguments "
                       "and result, plus ``cg_threads`` and ``cg_fallback`` (see pycalphad_cg.Workspace).\n\n"
                       + (_pycalphad_equilibrium.__doc__ or ''))
