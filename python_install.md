# Installing the python wrappers

First you must have the python module of ssht installed (`pyssht`). This can be found at http://astro-informatics.github.io/ssht/.

For some of the functions healpy must be installed (https://healpy.readthedocs.io/en/latest/). This is not necessary for all the code so can be skipped if you do not want to use the HEALPix sampling scheme.

To install (if you have healpy) run:

`python setup.py build_ext --inplace`

To install (if you do not have healpy) run:

`python setup_no_healpy.py build_ext --inplace`
