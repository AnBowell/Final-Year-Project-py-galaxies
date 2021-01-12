from distutils.core import setup
#import setuptools
from Cython.Build import cythonize
from distutils.extension import Extension

sourcefiles = ['c_test.pyx', 'simple_c_test.c']

extensions = [Extension("c_test", sourcefiles)]

setup(
    ext_modules = cythonize(extensions)
)
