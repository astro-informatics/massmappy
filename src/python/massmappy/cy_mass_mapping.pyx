# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
import pyssht as ssht
cimport cy_mass_mapping as mm
from libc.math cimport log, exp, sqrt
import matplotlib.pyplot as plt

#----------------------------------------------------------------------------------------------------#

def lm2lm_hp(np.ndarray[double complex, ndim=1, mode="c"] f_lm not None, int L):

        cdef np.ndarray[complex, ndim=1] f_lm_hp
        cdef int el, em, index

        f_lm_hp = np.empty([L*(L+1)/2,], dtype=complex)
        
        for el from 0 <= el < L:
                        for em from 0 <= em <= el:
                                        index = mm.cy_healpy_lm2ind(el, em, L)
                                        f_lm_hp[index] = f_lm[ el * el + el + em ]

        return f_lm_hp

#----------------------------------------------------------------------------------------------------#

def lm_hp2lm(np.ndarray[double complex, ndim=1, mode="c"] flm_hp not None, int L):

        cdef np.ndarray[complex, ndim=1] f_lm
        cdef int el, em, index

        f_lm = np.empty([L*L,], dtype=complex)
        for el from 0 <= el < L:
                        f_lm[ el * el + el] = flm_hp[mm.cy_healpy_lm2ind(el, 0, L)]  # m=0 case
                        for em from 1 <= em <= el:
                                        index = mm.cy_healpy_lm2ind(el, em, L)
                                        f_lm[ el * el + el - em ] = pow(-1.0, -em) * (flm_hp[index]).conjugate()
                                        f_lm[ el * el + el + em ] = flm_hp[index]

        return f_lm

#----------------------------------------------------------------------------------------------------#

def healpy_lm2ind(int el, int em, int L):
        return em*(2*L-1-em)/2+el

def healpy_ind2lm(int ind, int L):

    cdef int n_elements = L, m=0, el

    while (True):
        if ind >= n_elements:
            ind -= n_elements
            n_elements -= <int>1
            m += 1
        else:
            el = ind + m
            return (el, m)



#--------------------------------------------------------------#

def generate_kappa_lm_hp(np.ndarray[double, ndim=1, mode="c"] Cl not None, int L, int seed=-1):
# generate converage harmonic coefs in HEALPix lm convention
    cdef np.ndarray[complex, ndim=1] k_lm
    cdef int el, em, index

    k_lm = np.empty((L*(L+1)/2,), dtype=complex)

    if seed > -1:
        np.random.seed(seed)

    k_lm[mm.cy_healpy_lm2ind(0, 0, L)] = 0.0
    k_lm[mm.cy_healpy_lm2ind(1, 0, L)] = 0.0
    k_lm[mm.cy_healpy_lm2ind(1, 1, L)] = 0.0

    for el in range(2,L):
        index = mm.cy_healpy_lm2ind(el, 0, L)
        k_lm[index] = np.random.randn()*sqrt(Cl[el])
        for em in range(1,el+1):
            index = mm.cy_healpy_lm2ind(el, em, L)
            k_lm[index] = (np.random.randn()+ 1j*np.random.randn())*sqrt(Cl[el])
    
    return k_lm
    
#--------------------------------------------------------------#

def generate_kappa_lm_mw(np.ndarray[double, ndim=1, mode="c"] Cl not None, int L, int seed=-1):
# generate converage harmonic coefs in SSHT lm convention
    cdef np.ndarray[complex, ndim=1] k_lm
    cdef int el, em, index, index2

    k_lm = np.empty((L*L,), dtype=complex)

    if seed > -1:
        np.random.seed(seed)

    k_lm[0] = 0.0; k_lm[1] = 0.0; k_lm[2] = 0.0; k_lm[3] = 0.0;

    for el in range(2,L):
        index = mm.cy_mw_elm2ind(el, 0)
        k_lm[index] = np.random.randn()*sqrt(Cl[el])
        for em in range(1,el+1):
            index  = mm.cy_mw_elm2ind(el, em)
            k_lm[index] = (np.random.randn()+ 1j*np.random.randn())*sqrt(Cl[el])
            index2 = mm.cy_mw_elm2ind(el,-em)
            k_lm[index2] = pow(-1.0, -em) * (k_lm[index] ).conjugate()
    return k_lm

#--------------------------------------------------------------#


def kappa_lm_to_gamma_lm_mw(np.ndarray[double complex, ndim=1, mode="c"] k_lm not None, int L):
# converting converace to shear in harmonic space in SSHT lm convention
    cdef np.ndarray[complex, ndim=1] gamma_lm
    cdef int el, em, index
    cdef float D_ell

    gamma_lm = np.empty(L*L, dtype=np.complex_)

    gamma_lm[0] = 0.0; gamma_lm[1] = 0.0; gamma_lm[2] = 0.0; gamma_lm[3] = 0.0
    
    for ell in range(2,L):
        D_ell = sqrt((<float>ell+2.0)*(<float>ell-1.0)/((<float>ell+1.0)*<float>ell))
        for em in range(-ell,ell+1):
            index = mm.cy_mw_elm2ind(ell, em)
            gamma_lm[index] = <complex>D_ell * k_lm[index]
            
    return gamma_lm

#--------------------------------------------------------------#


def kappa_lm_to_gamma_lm_hp(np.ndarray[double complex, ndim=1, mode="c"] k_lm not None, int L):
# converting converace to shear in harmonic space in healpy lm convention
    cdef np.ndarray[complex, ndim=1] gamma_E_lm, gamma_B_lm
    cdef int el, em, index
    cdef float D_ell

    gamma_E_lm = np.empty((L*(L+1)/2,), dtype=complex)
    gamma_B_lm = np.zeros((L*(L+1)/2,), dtype=complex)

    gamma_E_lm[mm.cy_healpy_lm2ind(0, 0, L)] = 0.0; 
    gamma_E_lm[mm.cy_healpy_lm2ind(1, 0, L)] = 0.0; 
    gamma_E_lm[mm.cy_healpy_lm2ind(1, 1, L)] = 0.0;

    for ell in range(2,L):
        D_ell = sqrt((<float>ell+2.0)*(<float>ell-1.0)/((<float>ell+1.0)*<float>ell))
        for em in range(0,ell+1):
            index = mm.cy_healpy_lm2ind(ell, em, L)
            gamma_E_lm[index] = <complex>D_ell * k_lm[index]
            
    return gamma_E_lm, gamma_B_lm

#--------------------------------------------------------------#

def gamma_lm_to_kappa_lm_mw(np.ndarray[double complex, ndim=1, mode="c"] gamma_lm not None, int L, float sigma=-1):

    cdef np.ndarray[complex, ndim=1] kappa_lm
    cdef int el, em, index
    cdef float D_ell

    kappa_lm = np.empty(L*L, dtype=np.complex_)

    kappa_lm[0] = 0.0; kappa_lm[1] = 0.0; kappa_lm[2] = 0.0; kappa_lm[3] = 0.0
    
    for ell in range(2,L):
        D_ell = sqrt((<float>ell+2.0)*(<float>ell-1.0)/((<float>ell+1.0)*<float>ell))
        if sigma > 0.0:
            guassian = exp(-<float>ell*<float>ell*sigma*sigma)
        for em in range(-ell,ell+1):
            index = mm.cy_mw_elm2ind(ell, em)
            kappa_lm[index] = gamma_lm[index]/<complex>D_ell
            if sigma > 0.0:
                kappa_lm[index] = guassian*kappa_lm[index]
            
    return kappa_lm

#--------------------------------------------------------------#

def gamma_lm_to_kappa_lm_hp(np.ndarray[double complex, ndim=1, mode="c"] gamma_E_lm not None, \
    np.ndarray[double complex, ndim=1, mode="c"] gamma_B_lm not None, int L, float sigma=-1):

    cdef np.ndarray[complex, ndim=1] kappa_E_lm, kappa_B_lm
    cdef int el, em, index
    cdef float D_ell, guassian

    kappa_E_lm = np.empty((L*(L+1)/2,), dtype=np.complex_)
    kappa_B_lm = np.empty((L*(L+1)/2,), dtype=np.complex_)

    kappa_E_lm[mm.cy_healpy_lm2ind(0, 0, L)] = 0.0; 
    kappa_E_lm[mm.cy_healpy_lm2ind(1, 0, L)] = 0.0; 
    kappa_E_lm[mm.cy_healpy_lm2ind(1, 1, L)] = 0.0;
    
    kappa_B_lm[mm.cy_healpy_lm2ind(0, 0, L)] = 0.0; 
    kappa_B_lm[mm.cy_healpy_lm2ind(1, 0, L)] = 0.0; 
    kappa_B_lm[mm.cy_healpy_lm2ind(1, 1, L)] = 0.0;
    
    for ell in range(2,L):
        D_ell = sqrt((<float>ell+2.0)*(<float>ell-1.0)/((<float>ell+1.0)*<float>ell))
        if sigma > 0.0:
            guassian = exp(-<float>ell*<float>ell*sigma*sigma)
        for em in range(0,ell+1):
            index = mm.cy_healpy_lm2ind(ell, em, L)
            kappa_E_lm[index] = gamma_E_lm[index]/<complex>D_ell
            kappa_B_lm[index] = gamma_B_lm[index]/<complex>D_ell
            if sigma > 0.0:
                kappa_E_lm[index] = kappa_E_lm[index]*guassian
                kappa_B_lm[index] = kappa_B_lm[index]*guassian

    return kappa_E_lm, kappa_B_lm

def reduced_shear_to_kappa_mw(np.ndarray[complex, ndim=2, mode="c"] gamma not None, int L, str Method="MW", float sigma=-1,\
    float tol_error=1E-10, bint Iterate=True, bint return_count=False):

    cdef np.ndarray[complex, ndim=1] gamma_lm, k_lm
    cdef np.ndarray[long, ndim=2] mask
    cdef np.ndarray[complex, ndim=2] gamma_smooth, gamma_dum, k_mw_1, k_mw_2
    cdef int i, j, n_theta, n_phi, count
    cdef bint rel_error=True

    n_theta, n_phi = ssht.sample_shape(L,Method=Method)

    gamma_lm = np.zeros((L*L), dtype=complex)
    kappa_lm = np.zeros((L*L), dtype=complex)

    mask         = np.full((n_theta,n_phi),1,dtype=int)
    gamma_smooth = np.zeros((n_theta,n_phi),dtype=complex)
    gamma_dum    = np.zeros((n_theta,n_phi),dtype=complex)
    k_mw_1       = np.zeros((n_theta,n_phi),dtype=complex)
    k_mw_2       = np.zeros((n_theta,n_phi),dtype=complex)

    for i in range(n_theta):
        for j in range(n_phi):
            if np.isnan(gamma[i,j]):
                mask[i,j] = 0
                gamma[i,j] = 0.0

    if sigma>0:
        gamma_lm     = ssht.forward(gamma, L, Method=Method, Spin=2)
        gamma_lm     = ssht.guassian_smoothing(gamma_lm, L, sigma_in=sigma)
        gamma_smooth = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)
    else:
        gamma_smooth = gamma

    k_mw_1 = cy_gamma_to_kappa_mw(gamma_smooth, gamma_lm, kappa_lm, k_mw_1, L, Method=Method, sigma=-1)

    if Iterate:
        count = 0
        while(rel_error):
            count += 1
            if count>500:
                raise RuntimeError('Max iterations reached')

            for i in range(n_theta):
                for j in range(n_phi):
                    gamma_dum[i,j] = gamma_smooth[i,j]*(1-k_mw_1[i,j].real)

            k_mw_2 = cy_gamma_to_kappa_mw(gamma_dum, gamma_lm, kappa_lm, k_mw_2, L, Method=Method, sigma=-1)

            rel_error = False
            for i in range(n_theta):
                for j in range(n_phi):
                    if np.abs(k_mw_2[i,j] -k_mw_1[i,j]) > tol_error:
                        rel_error = True

            for i in range(n_theta):
                for j in range(n_phi):
                    k_mw_1[i,j] = k_mw_2[i,j]


    for i in range(n_theta):
        for j in range(n_phi):
            if mask[i,j] == 0:
                gamma[i,j] = np.nan + 1j*np.nan
                k_mw_1[i,j] = np.nan + 1j*np.nan

    if return_count:
        return k_mw_1, count
    else:
        return k_mw_1


cdef cy_gamma_to_kappa_mw(np.ndarray[complex, ndim=2, mode="c"] gamma, np.ndarray[complex, ndim=1, mode="c"] gamma_lm,\
    np.ndarray[complex, ndim=1, mode="c"] kappa_lm, np.ndarray[complex, ndim=2, mode="c"] k_mw,\
    int L, str Method="MW", float sigma=-1):

    cdef int i, j, n_theta, n_phi

    n_theta, n_phi = ssht.sample_shape(L,Method=Method)

    gamma_lm = ssht.forward(gamma, L, Method=Method, Spin=2)

    k_lm = gamma_lm_to_kappa_lm_mw(gamma_lm, L, sigma=sigma)

    k_mw = ssht.inverse(k_lm, L, Method=Method)

    return k_mw

def gamma_to_kappa_mw(np.ndarray[complex, ndim=2, mode="c"] gamma not None, int L, str Method="MW", float sigma=-1):

    cdef np.ndarray[complex, ndim=1] gamma_lm, k_lm
    cdef np.ndarray[long, ndim=2] mask
    cdef np.ndarray[complex, ndim=2] k_mw
    cdef int i, j, n_theta, n_phi

    n_theta, n_phi = ssht.sample_shape(L,Method=Method)

    mask = np.full((n_theta,n_phi),1,dtype=int)

    for i in range(n_theta):
        for j in range(n_phi):
            if np.isnan(gamma[i,j]):
                mask[i,j] = 0
                gamma[i,j] = 0.0

    gamma_lm = ssht.forward(gamma, L, Method=Method, Spin=2)

    k_lm = gamma_lm_to_kappa_lm_mw(gamma_lm, L, sigma=sigma)

    k_mw = ssht.inverse(k_lm, L, Method=Method)

    for i in range(n_theta):
        for j in range(n_phi):
            if mask[i,j] == 0:
                gamma[i,j] = np.nan + 1j*np.nan
                k_mw[i,j] = np.nan + 1j*np.nan


    return k_mw

def reduced_shear_to_kappa_plane(np.ndarray[complex, ndim=2, mode="c"] gamma not None, float delta_theta, float delta_phi, \
    float sigma=-1,float tol_error=1E-10, bint Iterate=True, bint return_count=False):

    cdef np.ndarray[complex, ndim=2] gamma_kk, k_kk, k_mw_1, k_mw_2, gamma_dum, gamma_smooth
    cdef np.ndarray[np.float_t, ndim=2] mask
    cdef int i, j, N, M, count
    cdef bint rel_error = True

    N = gamma.shape[0]
    M = gamma.shape[1]

    mask         = np.zeros((N,M), dtype=float)
    k_mw_1       = np.zeros((N,M), dtype=complex)
    k_mw_2       = np.zeros((N,M), dtype=complex)     
    gamma_dum    = np.zeros((N,M), dtype=complex)     
    gamma_kk     = np.zeros((N,M), dtype=complex)     
    gamma_smooth = np.zeros((N,M), dtype=complex)     
    k_kk         = np.zeros((N,M), dtype=complex)     

    for i in range(N):
        for j in range(M):
            if np.isnan(gamma[i,j]):
                gamma[i,j] = 0.0 + 0.0j
                mask[i,j]  = np.nan

    if sigma>0:
        gamma_kk = np.fft.fft2(gamma,norm="ortho")
        gamma_kk = np.fft.fftshift(gamma_kk)
        for i in range(N):
            for j in range(M):
                l1 = (<float>i-<float>(N/2))/delta_theta; 
                l2 = (<float>N*(<float>j-<float>(M/2)))/(<float>M*delta_phi);
                gamma_kk[i,j] = gamma_kk[i,j]*exp(-((l1*l1*sigma*sigma/(N*N))+(l2*l2*sigma*sigma/(M*M)))*np.pi*np.pi)
        gamma_kk = np.fft.ifftshift(gamma_kk)
        gamma_smooth = np.fft.ifft2(gamma_kk,norm="ortho")
    else:
        gamma_smooth = gamma


    gamma_kk = np.fft.fft2(gamma_smooth,norm="ortho")

    gamma_kk = np.fft.fftshift(gamma_kk)

    for i in range(N):
        for j in range(M):
            l1 = (<float>i-<float>(N/2))/delta_theta; 
            l2 = (<float>N*(<float>j-<float>(M/2)))/(<float>M*delta_phi);
        
            if not(abs(l1) < 1E-6 and abs(l2) < 1E-6):
                D_ij = (l1*l1 - l2*l2 - 2*1j*l1*l2)/(l1*l1+l2*l2)
                k_kk[i,j] = gamma_kk[i,j]*D_ij
            else:
                k_kk[i,j] = 0.0

    k_kk = np.fft.ifftshift(k_kk)

    k_mw_1 = np.fft.ifft2(k_kk,norm="ortho")

    if Iterate:
        count = 0
        while(rel_error):
            count += 1
            if count>500:
                raise RuntimeError('Max iterations reached')

            for i in range(N):
                for j in range(M):
                    gamma_dum[i,j] = gamma_smooth[i,j]*(1-k_mw_1[i,j].real)

            gamma_kk = np.fft.fft2(gamma_dum,norm="ortho")

            gamma_kk = np.fft.fftshift(gamma_kk)

            for i in range(N):
                for j in range(M):
                    l1 = (<float>i-<float>(N/2))/delta_theta; 
                    l2 = (<float>N*(<float>j-<float>(M/2)))/(<float>M*delta_phi);
            
                    if not(abs(l1) < 1E-6 and abs(l2) < 1E-6):
                        D_ij = (l1*l1 - l2*l2 - 2*1j*l1*l2)/(l1*l1+l2*l2)
                        k_kk[i,j] = gamma_kk[i,j]*D_ij
                    else:
                        k_kk[i,j] = 0.0

            k_kk = np.fft.ifftshift(k_kk)

            k_mw_2 = np.fft.ifft2(k_kk,norm="ortho")

            rel_error = False
            for i in range(N):
                for j in range(M):
                    if np.abs(k_mw_2[i,j] -k_mw_1[i,j]) > tol_error:
                        rel_error = True

            for i in range(N):
                for j in range(M):
                    k_mw_1[i,j] = k_mw_2[i,j]


    for i in range(N):
        for j in range(M):
            if np.isnan(mask[i,j]):
                gamma[i,j] = np.nan + 1j*np.nan
                k_mw_1[i,j] = np.nan + 1j*np.nan

    if return_count:
        return k_mw_1, count
    else:
        return k_mw_1


def gamma_to_kappa_plane(np.ndarray[complex, ndim=2, mode="c"] gamma not None, float delta_theta, float delta_phi, float sigma=-1):

    cdef np.ndarray[complex, ndim=2] gamma_kk, k_kk, k_mw
    cdef np.ndarray[np.float_t, ndim=2] mask
    cdef int i, j, N, M
    cdef complex D_ij
    cdef float l1, l2

    N = gamma.shape[0]
    M = gamma.shape[1]

    mask = np.zeros((N,M), dtype=float)

    for i in range(N):
        for j in range(M):
            if np.isnan(gamma[i,j]):
                gamma[i,j] = 0.0 + 0.0j
                mask[i,j]  = np.nan

    gamma_kk = np.fft.fft2(gamma,norm="ortho")

    gamma_kk = np.fft.fftshift(gamma_kk)

    k_kk = np.empty((N,M), dtype=complex)

    for i in range(N):
        for j in range(M):
            l1 = (<float>i-<float>(N/2))/delta_theta; 
            l2 = (<float>N*(<float>j-<float>(M/2)))/(<float>M*delta_phi);
#            l1 = (<float>i-<float>(N/2)) 
#            l2 = (<float>N*(<float>j-<float>(M/2)))/(<float>M);
        
            if not(abs(l1) < 1E-6 and abs(l2) < 1E-6):
                D_ij = (l1*l1 - l2*l2 - 2*1j*l1*l2)/(l1*l1+l2*l2)
                k_kk[i,j] = gamma_kk[i,j]*D_ij
                if sigma > 0.0:
                    k_kk[i,j] = k_kk[i,j]*exp(-((l1*l1*sigma*sigma/(N*N))+(l2*l2*sigma*sigma/(M*M)))*np.pi*np.pi)
#                    print l1, sigma,  l1*l1*sigma*sigma/N
            else:
                k_kk[i,j] = 0.0

    k_kk = np.fft.ifftshift(k_kk)

    k_mw = np.fft.ifft2(k_kk,norm="ortho")

    for i in range(N):
        for j in range(M):
            if np.isnan(mask[i,j]):
                gamma[i,j] = np.nan + 1j*np.nan
                k_mw[i,j] = np.nan + 1j*np.nan

    return k_mw
