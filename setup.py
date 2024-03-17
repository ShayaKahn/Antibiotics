from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(['GLV_functions.pyx', 'Null_model_functions.pyx'], annotate=True),
    include_dirs=[numpy.get_include()]
)