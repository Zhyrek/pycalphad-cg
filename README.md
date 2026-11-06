# pycalphad-cg

[pycalphad](https://github.com/pycalphad/pycalphad) with a column-generation equilibrium solver, as a drop-in
replacement for the `pycalphad` namespace. It installs alongside pycalphad and needs no changes to it.

```python
import pycalphad_cg as pyc          # in place of: import pycalphad as pyc
from pycalphad_cg import Database, equilibrium, binplot, variables as v

dbf = Database('alfe.tdb')
conds = {v.T: (500, 2000, 10), v.X('AL'): (0.01, 0.99, 0.01), v.P: 101325, v.N: 1}
eq = equilibrium(dbf, ['AL', 'FE', 'VA'], list(dbf.phases), conds, output='HM')
ax = binplot(dbf, ['AL', 'FE', 'VA'], list(dbf.phases), conds)
```

Everything that computes equilibria uses column generation: `equilibrium`, `Workspace`, the mapping strategies
(`pycalphad_cg.mapping.BinaryStrategy` etc.), `binplot` and `ternplot`. Everything else is pycalphad's own object
(`pycalphad_cg.Database is pycalphad.Database`), and so are pycalphad's modules: `pycalphad_cg.X` is `pycalphad.X`
for every module this package does not replace, so `import pycalphad_cg.variables as v` or
`from pycalphad_cg.plot.utils import phase_legend` work too. (The import name has an underscore: Python names
cannot contain hyphens.)

## Install

From PyPI, once released (see [RELEASING.md](RELEASING.md)):

```bash
pip install pycalphad-cg
```

From this repository (needs a C++ compiler):

```bash
pip install git+https://github.com/Zhyrek/pycalphad-cg
```

It works with pycalphad 0.11.2 (the current release) and pycalphad's develop branch. The compiled core uses no C
declarations from pycalphad (it reads pycalphad's compiled energy functions through its Python API, the same way
pycalphad itself does), so one build works with any pycalphad version that keeps that API.

## How column generation finds an equilibrium

pycalphad samples every phase on a grid, takes the lower convex hull of the samples as starting points, and refines
them by Newton's method. Column generation needs no grid:

- **Linear program.** The equilibrium over a set of points (a phase at given site fractions) is a linear program.
  Its dual values are the chemical potentials: a hyperplane.
- **Pricing.** Each round minimizes, for every phase, its energy relative to the current hyperplane (the
  tangent-plane distance), by local minimization from several starts. Every minimum below the hyperplane becomes a
  new point, and the linear program is solved again.
- **Polish.** When the phases are identified, Newton's method finishes the equilibrium on its defining equations.
- **Certification.** A final pricing round at the polished chemical potentials confirms that no phase lies below
  the hyperplane.

Because the minimization works directly against the hyperplane, column generation finds states that a fixed grid
can miss (for example ordered states that occupy a small region of site-fraction space). Conditions that share
state variables are solved in sequence, each starting from the previous one's points, on worker threads that run
without the GIL.

## Options

```python
pyc.options.threads = 4       # worker threads (default: PYCALPHAD_CG_THREADS, else all cores)
pyc.options.fallback = True   # unsupported problems: pycalphad's solver with a warning (False: raise)

eq = pyc.equilibrium(dbf, comps, phases, conds, cg_threads=1)         # per call
wks = pyc.Workspace(dbf, comps, phases, conds, cg_fallback=False)     # per Workspace
```

`pyc.cg_mode()` applies column generation to pycalphad's own functions within a `with` block (in that thread
only):

```python
import pycalphad
with pyc.cg_mode():
    eq = pycalphad.equilibrium(dbf, comps, phases, conds)
```

## How it plugs into pycalphad

- **`Workspace`** subclasses pycalphad's and overrides `recompute()` only, which returns the equilibrium dataset in
  pycalphad's layout. Everything that reads it (`wks.get(...)`, computed properties, Jansson derivatives, reference
  states) is pycalphad's code, unchanged.
- **`equilibrium` and mapping.** pycalphad has no hook for an outside solver here: `pycalphad.equilibrium` and the
  mapping strategies' starting points create a `Workspace`, and the mapping steps create pycalphad's `Solver` and
  search a grid for new phases, each by a module-level name. While a thread is inside `cg_mode()` (which
  `pyc.equilibrium`, `pyc.binplot`, `pyc.ternplot` and the methods of the strategies in `pyc.mapping` enter), those
  four names resolve to the column-generation versions for that thread. Other threads keep getting pycalphad's own,
  and the names are restored when no thread is in the mode.
- **Plugin import.** pycalphad imports every installed `pycalphad_*` package when it is imported. pycalphad-cg
  changes nothing at import, so plain `pycalphad` behaves as before when pycalphad-cg is installed.

The pycalphad internals this relies on are listed under [Compatibility](#compatibility).

## Validation

**pycalphad's own test suite with column generation.** `tests/upstream_suite_plugin.py` runs pycalphad's tests
with every test inside `cg_mode()` and `pycalphad.Workspace` replaced:

| pycalphad | plain | with column generation | column-generation solves | fallbacks |
|---|---|---|---|---|
| 0.11.2 | 309 passed, 2 skipped, 1 xfailed | 309 passed, 2 skipped, 1 xfailed | 173 | 5 |
| develop (0cbc613c0) | DEV_PLAIN | DEV_CG | DEV_SOLVES | DEV_FALLBACKS |

The fallbacks are the tests that pass their own solver object (see Limits).

**This package's tests** (`pytest tests`) compare with pycalphad on binaries and ternaries from pycalphad's test
databases:

- column generation's Gibbs energy is never above pycalphad's by more than 0.01 J/mol, and where the two agree the
  chemical potentials do too;
- `Workspace.get` properties, fixed chemical potentials and mass fractions, and the `equilibrium` output match;
- input errors match; fallback works; 1 and 4 threads give the same result;
- mapping goes through column generation, and plain pycalphad in the same process (also in another thread at the
  same time) does not;
- the drop-in namespace, module aliases and both import orders.

## Performance

PERFORMANCE_TABLE

## Limits

- **Fallback to pycalphad's solver** (with a warning; `fallback=False` raises `CGCapabilityError` instead):
  - `calc_opts` other than `pdens`;
  - custom solver objects (pycalphad's grid-and-hull starting points would be skipped, so that solver would never
    see them) and custom phase-record factories;
  - a phase without a temperature state variable, or with non-linear internal constraints;
  - conditions other than T, P, N, mole and mass fractions, chemical potentials and linear combinations of mole
    fractions;
  - phase-local conditions other than site fractions and mole fractions, or a phase-local mole fraction in a phase
    whose formula moles depend on its site fractions (ionic liquids).
- **Properties that solve their own sub-problems** (`DormantPhase`, `T0`, driving forces through
  `pycalphad.property_framework`) use pycalphad's Newton solver on the composition sets the equilibrium returned,
  as in pycalphad.
- **Deep module paths** give pycalphad's objects: `pycalphad_cg.mapping.BinaryStrategy` and
  `pycalphad_cg.mapping.strategy.BinaryStrategy` use column generation, but
  `pycalphad_cg.mapping.strategy.binary_strategy.BinaryStrategy` is pycalphad's own class.
- **One calculation at a time per process.** The compiled phases are cached across calls and shared, so concurrent
  calls from several threads run one after the other (each uses its own worker threads).

## Compatibility

Nothing here needs changes to pycalphad, but these internals can change in a pycalphad release (CI runs against
pycalphad's develop branch weekly to notice):

- `Workspace.recompute()` as the override point, and the `LightDataset` layout of `wks.eq` (variables GM, MU, NP, X,
  Y, Phase; `vertex` of length C+1; Y padded to the largest phase; empty conditions left NaN);
- `PhaseRecordFactory.get_phase_property`, `get_phase_formula_moles_element`, `update_parameters`, `param_values`,
  `param_symbols`, `molar_masses`, `nonvacant_elements`, `state_variables`, `models`, and SymEngine's `as_ctypes()`
  on the compiled functions;
- `PhaseRecord`'s Python attributes and its internal-constraint and formula-mole methods;
- the module-level names `pycalphad.core.equilibrium.Workspace`, `pycalphad.mapping.starting_points.Workspace`,
  `pycalphad.mapping.zpf_equilibrium.Solver` and `pycalphad.mapping.zpf_equilibrium._find_global_min_cs` (if one is
  missing, pycalphad-cg warns and that path keeps pycalphad's solver);
- the input checks of `calculate` and `starting_point`, re-implemented with the same errors;
- `CompositionSet` and `SolverBase`/`SolverResult` for the mapping solver.

## Development

```bash
pip install --no-build-isolation -e .     # needs Cython and SciPy
pytest tests
# pycalphad's own suite with column generation:
cd /tmp && PYTHONPATH=/path/to/pycalphad-cg/tests pytest --pyargs pycalphad.tests -p upstream_suite_plugin
```

Written with the help of an LLM (Claude); every change was reviewed by the author.

## License

MIT, like pycalphad.
