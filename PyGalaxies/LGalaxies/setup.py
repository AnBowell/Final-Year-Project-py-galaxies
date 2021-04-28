from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension



sourcefiles = ['L_Galaxies.pyx', 'cool_func.c','c_model_params.c',
               'model_infall.c','model_star_formation_and_feedback.c']


extensions = [Extension("L_Galaxies", sourcefiles)]


setup(
    ext_modules = cythonize(extensions)

)


