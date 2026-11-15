"""census.py with every map inside pycalphad_cg.cg_mode() (pycalphad's mapping, column generation's solves)."""
import runpy, sys, os
import pycalphad_cg
sys.argv = [os.path.join(os.path.dirname(__file__), 'census.py')] + sys.argv[1:]
with pycalphad_cg.cg_mode():
    runpy.run_path(sys.argv[0], run_name='__main__')
