from distutils.core import setup
from Cython.Build import cythonize
 
setup(
	name = 'MHM_elastodynamics app',
	ext_modules = cythonize("MHM_elastodynamics.pyx"),
)
