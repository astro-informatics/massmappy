import sys
import os
import shutil

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

# clean previous build
for root, dirs, files in os.walk("./src/python/", topdown=False):
    for name in dirs:
        if (name == "build"):
            shutil.rmtree(name)


include_dirs = [
    numpy.get_include(),
   os.environ['S2LET']+"/include/",
   os.environ['SO3']+"/include/c/",
   os.environ['SSHT']+"/include/c/",
    ]

extra_link_args=[
   "-L"+os.environ['S2LET']+"/lib/",
   "-L"+os.environ['SO3']+"/lib/c/",
   "-L"+os.environ['SSHT']+"/lib/c/",
   "-L"+os.environ['FFTW']+"/lib",

]

setup(
    classifiers=['Programming Language :: Python :: 2.7'],
    name = "pyssht",
    version = "2.0",
    prefix='.',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize([Extension(
        "src/python/cy_mass_mapping",
        package_dir=['src'],
        sources=["src/python/cy_mass_mapping.pyx"],
        include_dirs=include_dirs,
        libraries=[],#["ssht", "fftw3"],
        extra_link_args=extra_link_args,
        extra_compile_args=[]
    ),
    Extension("src/python/cy_healpy_mass_mapping",
            package_dir=['src'],
            sources=["src/python/cy_healpy_mass_mapping.pyx"],
            include_dirs=include_dirs,
            libraries=[],
            extra_link_args=extra_link_args,
            extra_compile_args=[]
    ),
    Extension("src/python/cy_lsst_challenge_1_utils",
            package_dir=['src'],
            sources=["src/python/cy_lsst_challenge_1_utils.pyx"],
            include_dirs=include_dirs,
            libraries=[],
            extra_link_args=extra_link_args,
            extra_compile_args=[]
    ),
    Extension("src/python/cy_eb_sep",
            package_dir=['src'],
            sources=["src/python/cy_eb_sep.pyx"],
            include_dirs=include_dirs,
            libraries=["s2let", "so3", "fftw3", "ssht"],
            extra_link_args=extra_link_args,
            extra_compile_args=[]
    ),
    Extension("src/python/cy_DES_utils",
            package_dir=['src'],
            sources=["src/python/cy_DES_utils.pyx"],
            include_dirs=include_dirs,
            libraries=[],
            extra_link_args=extra_link_args,
            extra_compile_args=[])
    ])
)