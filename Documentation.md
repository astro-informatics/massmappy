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



