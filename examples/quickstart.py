"""Al-Fe with pycalphad_cg, used exactly like pycalphad: only the import changes."""
from importlib.resources import files

import matplotlib
import numpy as np

import pycalphad_cg as pyc   # in place of: import pycalphad as pyc
from pycalphad_cg import variables as v

matplotlib.use('Agg')

dbf = pyc.Database(str(files('pycalphad.tests.databases') / 'alfe.tdb'))
comps = ['AL', 'FE', 'VA']
phases = list(dbf.phases)
conds = {v.T: (600, 1801, 100), v.X('AL'): (0.05, 0.96, 0.05), v.P: 101325, v.N: 1}

# equilibrium: same arguments and result as pycalphad.equilibrium
eq = pyc.equilibrium(dbf, comps, phases, conds, output='HM')
print(eq)

# Workspace: the equilibria by column generation, everything after them (properties, derivatives) is pycalphad's
wks = pyc.Workspace(dbf, comps, phases, {**conds, v.T: 1000})
x_al, cp = wks.get('X(AL)', 'HM.T')
print('heat capacity at 1000 K:')
for x, c in zip(np.ravel(x_al)[::3], np.ravel(cp)[::3]):
    print(f'  x(Al) = {x:.2f}: {c:.2f} J/mol/K')

# phase diagram: pycalphad's mapping, with column generation for its equilibria and steps
ax = pyc.binplot(dbf, comps, phases, {v.T: (500, 2000, 10), v.X('AL'): (0, 1, 0.01), v.P: 101325, v.N: 1})
ax.figure.savefig('alfe.png', dpi=100)
print('phase diagram written to alfe.png')
