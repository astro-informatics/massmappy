# massmappy documentation

Here is a description of each function the inputs and outputs. If you have any quations then please email me at christophergrwallis@gmail.com and I'll help you out (and probably update this doc to make it clearer)

### lm2lm_hp

```lm2lm_hp(np.ndarray[double complex, ndim=1, mode="c"] f_lm not None, int L)```

Function to convert between SSHT and HEALPix definitions of how harmonic coefficiants are stored.

Inputs:
* `f_lm` array containing harmonic coefficiants using SSHT convention
* `L` Harmonic bandlimit

Outputs:
* Array containing array in HEALPix convention


### lm_hp2lm

```lm_hp2lm(np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, int L)```

Function to convert between HEALPix and SSHT definitions of how harmonic coefficiants are stored.

Inputs:
* `f_lm` array containing harmonic coefficiants using HEALPix convention
* `L` Harmonic bandlimit

Outputs:
* Array containing array in SSHT convention

### healpy_lm2ind

```healpy_lm2ind(int el, int em, int L):```

Function to convert from harmonic coefficiant terms ell and em to the index of a HEALPix array

Inputs:
* `el` harmonic identifier ell
* `em` harmonic identifier em
* `L` the harmonic band limit

Output:
* Integer of the index value

### healpy_ind2lm

```healpy_ind2lm(int ind, int L):```

Function to convert from the index of a HEALPix array to harmonic coefficiant terms ell and em 

Inputs:
* `ind` the index of the array
* `L` the harmonic band limit

Output:
* tuple containing:
    * `el` harmonic identifier ell
    * `em` harmonic identifier em


### generate_kappa_lm_hp

```generate_kappa_lm_hp(np.ndarray[double, ndim=1, mode="c"] Cl not None, int L, int seed=-1)```

Generate convergance harmonic coefs in HEALPix lm convention (assuming a Guassian realisation)

Inputs:
* `Cl` Array containing the E-mode power spectrum
* `L` the harmonic band limit
* `seed` seed for random number generator, set to -1 for different seed every time. (default=-1)
    
Output:
* Array containing harmonic coefficiants in HEALPix convention

### generate_kappa_lm_mw

```generate_kappa_lm_mw(np.ndarray[double, ndim=1, mode="c"] Cl not None, int L, int seed=-1):```

Generate converage harmonic coefs in SSHT lm convention


* `Cl` Array containing the E-mode power spectrum
* `L` the harmonic band limit
* `seed` seed for random number generator, set to -1 for different seed every time. (default=-1)
    
Output:
* Array containing convergence harmonic coefficiants in SSHT convention

### kappa_lm_to_gamma_lm_mw

```kappa_lm_to_gamma_lm_mw(np.ndarray[double complex, ndim=1, mode="c"] k_lm not None, int L):```

Converting converace to shear in harmonic space in SSHT lm convention

Inputs:
* `k_lm` convergence  harmonic coefficiants in SSHT lm convention
* `L` harmonic bandlimit

Output:
* Array containing shear harmonic coefficiants in SSHT convention

### kappa_lm_to_gamma_lm_hp

```kappa_lm_to_gamma_lm_hp(np.ndarray[double complex, ndim=1, mode="c"] k_lm not None, int L):```

Converting convergence to shear in harmonic space in HEALPix lm convention
 
Inputs:
* `k_lm` convergence  harmonic coefficiants in HEALPix lm convention
* `L` harmonic bandlimit

Output:
* tuple containing:
    * Array containing shear E-mode harmonic coefficiants in HEALPix convention
    * Array containing shear B-mode harmonic coefficiants in HEALPix convention

### gamma_lm_to_kappa_lm_mw

```gamma_lm_to_kappa_lm_mw(np.ndarray[double complex, ndim=1, mode="c"] gamma_lm not None, int L, float sigma=-1):```

Converting shear to convergence in harmonic space in SSHT lm convention

Inputs:
* `gamma_lm` shear harmonic coefficiants in SSHT lm convention
* `L` harmonic bandlimit
* `sigma` standard deviation of a smoothing kernal to be applied in radians, if -1 then no smoothing is applied (default=-1)

Output:
* Array containing convergence harmonic coefficiants in SSHT convention

### gamma_lm_to_kappa_lm_hp

```gamma_lm_to_kappa_lm_hp(np.ndarray[double complex, ndim=1, mode="c"] gamma_E_lm not None, \
    np.ndarray[double complex, ndim=1, mode="c"] gamma_B_lm not None, int L, float sigma=-1):```

Converting shear to convergence in harmonic space in HEALPix lm convention

Inputs:
* `gamma_E_lm` shear E-mode harmonic coefficiants in SSHT lm convention
* `gamma_B_lm` shear B-mode harmonic coefficiants in SSHT lm convention
* `L` harmonic bandlimit
* `sigma` standard deviation of a smoothing kernal to be applied in radians, if -1 then no smoothing is applied (default=-1)

Output:
* Array containing convergence harmonic coefficiants in SSHT convention

### reduced_shear_to_kappa_mw

```reduced_shear_to_kappa_mw(np.ndarray[complex, ndim=2, mode="c"] gamma not None, int L, str Method="MW", float sigma=-1,\
    float tol_error=1E-10, bint Iterate=True, bint return_count=False):```

Converting reduced shear to convergence in real space in SSHT convention using an iterative scheme.

Inputs:
* `gamma` reduced shear in SSHT convention
* `L` harmonic bandlimit
* `Method` SSHT sampling used (see SSHT documentation) (default='MW')
* `sigma` standard deviation of a smoothing kernal to be applied before iterating in radians, if -1 then no smoothing is applied (default=-1)
* `tol_error` difference between two iterative steps after which converence is considered to have been reached (default=1E-10)
* `Iterate` bolian to say if an iterations should be performed if `False` the input is treated as shear (default=`True`)
* `return_count` if true the number of iterations required is returned (default=`False`)

Output:
If `return_count=False`
* Array containing convergence harmonic map in SSHT convention
If `return_count=True`
* A tuple containing:
    * Array containing convergence harmonic map in SSHT convention
    * Integer containing the number of iterations


### gamma_to_kappa_mw

```gamma_to_kappa_mw(np.ndarray[complex, ndim=2, mode="c"] gamma not None, int L, str Method="MW", float sigma=-1)```

Converting shear to convergence in real space in SSHT convention.

Inputs:
* `gamma` reduced shear in SSHT convention
* `L` harmonic bandlimit
* `Method` SSHT sampling used (see SSHT documentation) (default='MW')
* `sigma` standard deviation of a smoothing kernal to be applied before iterating in radians, if -1 then no smoothing is applied (default=-1)

Output:
* Array containing convergence harmonic map in SSHT convention