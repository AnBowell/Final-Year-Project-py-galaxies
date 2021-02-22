from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension



sourcefiles = ['L_Galaxies.pyx', 'cool_func.c','c_model_params.c',
               'model_infall.c']


extensions = [Extension("L_Galaxies", sourcefiles)]


setup(
    ext_modules = cythonize(extensions)

)


