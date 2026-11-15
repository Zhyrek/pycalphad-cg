# Mapping census

Maps a fixed set of 12 binary phase diagrams and 5 ternary isothermal sections with pycalphad's mapping and records
every ZPF line and node, to compare mapping results and times between pycalphad and pycalphad-cg (or two versions),
and checks every map against point equilibria.

```bash
cd /tmp   # (a neutral directory)
PYTHONHASHSEED=0 python /path/to/tools/census/census.py plain.json      # pycalphad's own solver
PYTHONHASHSEED=0 python /path/to/tools/census/census_cg.py cg.json      # the same maps inside pycalphad_cg.cg_mode()
python /path/to/tools/census/compare.py plain=plain.json cg=cg.json     # lines, failures, duplicates; regions one lacks
python /path/to/tools/census/groundtruth.py plain=plain.json cg=cg.json # where each map shows the wrong phases
```

Mapping is deterministic only with a fixed `PYTHONHASHSEED`. Databases come from pycalphad's test databases; the
Ba-Ca map needs `BaCa-86Alc.tdb`, which is not included (pass its directory with `--dbdir`; without it that map is
skipped).

Phases are compared by state (`states.py`): an ordered phase in its disordered state is the disordered phase, and a
phase that is an exact copy of an ordered phase's disordered part (FCC_A1 beside FCC_4SL in the Co-V and Au-Cu
databases) is that state too. So two lines of one boundary under two names count as duplicates.

`groundtruth.py` probes each diagram on a grid (binaries: 40 temperatures x 100 compositions; ternary sections: a
0.02 composition grid), and at the middle of every tie-line one map draws and another lacks. At each probe it takes
the phases each map shows there (three-phase triangle, two-phase strip between consecutive tie-lines, or the single
phase met first by rays to the phase boundaries), computes the equilibrium restricted to those phases and the
unrestricted point equilibrium (pycalphad with pdens 1000, 500 for ternaries, and pycalphad-cg), and takes the
lowest Gibbs energy as the truth. A map is wrong at a probe if its phases' equilibrium is higher than the truth by
more than `--tol` J/mol (default 1) or has other phases than the map shows (an ordered phase and its disordered
state count as one: maps do not show ordering of second order). Point equilibria above the truth are reported per
solver. A full run takes about 30 minutes; `--only` takes a comma-separated list of map names.

Single probes next to a boundary are usually the boundary's position within the probe spacing; regions of several
probes, or energy errors of hundreds of J/mol, are mapping errors.
