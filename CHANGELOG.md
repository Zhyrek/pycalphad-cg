# Changelog

## 0.1.0 (unreleased)

First release.

- Column-generation equilibria for pycalphad, as a drop-in for the `pycalphad` namespace (`import pycalphad_cg as pyc`):
  `equilibrium`, `Workspace`, the mapping strategies (`TielineStrategy` and the others), `binplot` and `ternplot`.
  Everything else, and every module not replaced, is pycalphad's own.
- Works with stock pycalphad 0.11.2, 0.12.0 and its develop branch, without changes to pycalphad: the compiled core
  uses no C declarations from pycalphad.
- Problems column generation does not handle fall back to pycalphad's solver (`options.fallback`).
- pycalphad's own test suite passes with column generation in place (see the README for the two tests that differ).
