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
(`pycalphad_cg.mapping.TielineStrategy` etc.), `binplot` and `ternplot`. Everything else is pycalphad's own object
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

It works with pycalphad 0.11.2 and 0.12.0 and with pycalphad's develop branch. The compiled core uses no C
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
| 0.12.0 | 414 passed, 2 skipped, 2 xfailed, 1 xpassed | 414 passed, 2 skipped, 2 xfailed, 1 xpassed | 260 | 6 |
| develop (a7ea88fac) | 414 passed, 2 skipped, 2 xfailed, 1 xpassed | 414 passed, 2 skipped, 2 xfailed, 1 xpassed | 260 | 6 |

The fallbacks are the tests that pass their own or a configured solver object (see Limits). (Apple M4, Python 3.12, SymEngine 0.14.1.)

Two of pycalphad's 0.12 mapping tests come out differently, which the plugin records:

- `test_zpf_line_from_saddle_point_configuration_reaches_invariant` (Ni-Pt), marked by pycalphad as an expected
  failure of its solver (pycalphad#746: L1<sub>2</sub> left on a saddle point), passes with column generation and is
  run as an ordinary test.
- `test_automatic_starting_points_for_every_merged_field` (Au-Cu) assumes that a step map at 300 K steps over the
  narrow Au<sub>3</sub>Cu field; column generation resolves that field instead. It is marked as an expected failure,
  and what it checks (a starting point in each direction for every field, the Au-rich field mapped) is a test of
  this package.

pycalphad#762 (Al-Ni with a four-sublattice FCC, `AlNi-22Yan.tdb` attached to the issue): both of its reference tests,
where pycalphad misses L1<sub>2</sub> below the invariants, pass with `pycalphad_cg.Workspace` in place of
`Workspace`.

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

`examples/timing.py` on an Apple M4 (10 cores): wall time per call, warm (the first call, which compiles the
phases' energy functions, is not timed), including building the Workspace (its models and phase records), which
is the same for both. The last column is the largest amount by which column generation's Gibbs energy is below
pycalphad's; it is never above pycalphad's by more than 0.0014 J/mol.

| System | Conditions | pycalphad develop | cg, 1 thread | cg, 10 threads | pycalphad 0.11.2 | cg lower by, at most (J/mol) |
|---|---|---|---|---|---|---|
| Al-Ni (alni_dupin_2001), 1200 K | 500 | 0.17 s | 0.07 s | 0.04 s | 0.28 s | 7e-5 |
| Al-Fe (Al-Fe_sundman2009), 900 K | 500 | 1.27 s | 0.43 s | 0.16 s | 1.55 s | 39 |
| Al-Fe (alfe), 500–1900 K | 15 × 25 | 0.12 s | 0.05 s | 0.03 s | 0.19 s | 22 |
| Al-Cr-Ni (alcrni), 1273 K | 25 × 25 | 0.40 s | 0.15 s | 0.06 s | 0.62 s | 150 |
| COST507 Al-Mg-Si-Cu-Zn, 800 K | 100 | 0.69 s | 0.44 s | 0.44 s | 1.57 s | 4e-6 |
| COST507, 10 elements, 800 K | 100 | 6.81 s | 2.62 s | 2.67 s | 21.0 s | 5e-4 |

(cg's times are the same on both pycalphad versions; the 0.11.2 column is pycalphad's own solver there.)

- **Where column generation is lower,** pycalphad's default grid (`pdens` 60) misses a state, and the state is
  real: re-evaluating column generation's phases with `pycalphad.calculate` at its site fractions gives its
  energy. Examples: Al-Fe at 700 K, x(Al) = 0.14, the ordered state of B2_BCC (22 J/mol lower); Al-Cr-Ni at
  1273 K, x(Al) = 0.14, x(Cr) = 0.18, L1<sub>2</sub> + L1<sub>2</sub> where pycalphad finds one L1<sub>2</sub> (152 J/mol
  lower). On a 0.01 grid of Al-Cr-Ni at 1273 K, column generation is lower than pycalphad even at `pdens`
  1000 at 1224 of 5050 compositions, by up to 512 J/mol (B2 + liquid at x(Al) = 0.665 on the Al-Cr edge, where
  pycalphad finds B2 alone), and never higher.
- **Threads** pay off from a few hundred conditions; at 100 conditions the fixed cost of building the Workspace
  dominates.
- **Mapping** with column generation took 2.2 s against pycalphad's 2.9 s for an Al-Ni binary (1000–2000 K), and
  2.3 s against 3.6 s for the Al-Cr-Ni isotherm at 1273 K.

## Limits

- **Fallback to pycalphad's solver** (with a warning; `fallback=False` raises `CGCapabilityError` instead):
  - `calc_opts` other than `pdens`;
  - custom solver objects, or pycalphad's `Solver` with non-default settings (pycalphad's grid-and-hull starting
    points would be skipped, so that solver would never see them), and custom phase-record factories;
  - a phase without a temperature state variable, or with non-linear internal constraints;
  - conditions other than T, P, N, mole and mass fractions, chemical potentials and linear combinations of mole
    fractions;
  - phase-local conditions other than site fractions and mole fractions, or a phase-local mole fraction in a phase
    whose formula moles depend on its site fractions (ionic liquids).
- **`fallback=False` in mapping:** pycalphad's mapping catches every error from its solver and goes on, so the
  `CGCapabilityError` is raised again when the outermost `cg_mode` block (or `binplot`, a strategy method, ...)
  ends.
- **Two descriptions of one phase** (e.g. an ordered model in its disordered state and the disordered model) are
  recognized in mapping by equal composition and equal energy (within 5 J/mol) also 100 K away. Two different
  phases at the same composition whose entropies differ by less than about 0.05 J/mol/K would be taken for one.
- **Properties that solve their own sub-problems** (`DormantPhase`, `T0`, driving forces through
  `pycalphad.property_framework`) use pycalphad's Newton solver on the composition sets the equilibrium returned,
  as in pycalphad.
- **Deep module paths** give pycalphad's objects: `pycalphad_cg.mapping.TielineStrategy` and
  `pycalphad_cg.mapping.strategy.TielineStrategy` use column generation, but
  `pycalphad_cg.mapping.strategy.tieline_strategy.TielineStrategy` is pycalphad's own class.
- **One calculation at a time per process.** The compiled phases are cached across calls and shared, so concurrent
  calls from several threads run one after the other (each uses its own worker threads).

## Where results differ from pycalphad

- Gibbs energies lower than pycalphad's where its grid misses a state (see Performance).
- Ratio conditions such as `v.X('NI') / v.X('ZN')` in systems with more components than the ratio names: pycalphad's
  own solver takes the denominator's column by its position among the condition's symbols rather than by element,
  so its result does not meet the condition in some systems; column generation meets it.

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
