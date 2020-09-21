from distutils.core import Extension,setup
from Cython.Build import cythonize

extensions = [
	Extension("potential", ["potential.pyx"],
		define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
		),
]

setup(
	ext_modules = cythonize(extensions, 
		compiler_directives={'language_level' : "3"})
	)
