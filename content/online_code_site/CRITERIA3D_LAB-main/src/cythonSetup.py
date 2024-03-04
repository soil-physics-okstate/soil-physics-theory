# compile with:
# python cythonSetup.py build_ext --inplace

from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Criteria 3D solver',
    ext_modules=cythonize("solverC.pyx"),
)
