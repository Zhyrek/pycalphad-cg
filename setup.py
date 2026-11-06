from Cython.Build import cythonize
from setuptools import Extension, setup

# The core uses no C declarations of pycalphad or NumPy (only SciPy's LAPACK, through scipy.linalg.cython_lapack),
# so one build works with any pycalphad version.
extensions = [Extension('pycalphad_cg._core', sources=['pycalphad_cg/_core.pyx'], language='c++')]
setup(ext_modules=cythonize(extensions, compiler_directives={'language_level': 3}))
