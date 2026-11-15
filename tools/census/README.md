# Mapping census

Maps a fixed set of 12 binary phase diagrams and 5 ternary isothermal sections with pycalphad's mapping and records
every ZPF line and node, to compare mapping results and times between pycalphad and pycalphad-cg (or two versions).

```bash
cd /tmp   # (a neutral directory)
PYTHONHASHSEED=0 python /path/to/tools/census/census.py plain.json      # pycalphad's own solver
PYTHONHASHSEED=0 python /path/to/tools/census/census_cg.py cg.json      # the same maps inside pycalphad_cg.cg_mode()
python /path/to/tools/census/compare.py plain=plain.json cg=cg.json     # lines, failures, duplicates; regions one lacks
python /path/to/tools/census/truth.py plain.json cg.json plain cg       # disputed regions against point equilibria
```

Mapping is deterministic only with a fixed `PYTHONHASHSEED`. Databases come from pycalphad's test databases; the
Ba-Ca map needs `BaCa-86Alc.tdb`, which is not included (pass its directory with `--dbdir`; without it that map is
skipped). `truth.py` computes, at the midpoint of each disputed ZPF point, the lower-energy equilibrium of pycalphad
(`pdens` 1000) and pycalphad-cg, and reports whether the disputed phase pair is stable there.
