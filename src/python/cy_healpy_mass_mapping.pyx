# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
import pyssht as ssht
import healpy as hp
import cy_mass_mapping as mw_mm
cimport cy_mass_mapping as mw_mm


def gamma_to_kappa_hp_boris(np.ndarray[double, ndim=1, mode="c"] gamma_real not None, \
    np.ndarray[double, ndim=1, mode="c"] gamma_imag not None, int L, int Nside):

    cdef int lmax=L-1

    maps_hp = [np.zeros(hp.nside2npix(Nside)),gamma_real, gamma_imag]
    [dummy, gamma_E_lm_rec, gamma_B_lm_rec] = hp.map2alm(maps_hp, lmax=lmax)
    E_lm = gamma_E_lm_rec
    ellfac = np.zeros(gamma_E_lm_rec.shape)
    for el in range(L):
        for em in range(L):
           ellfac[mw_mm.healpy_lm2ind(el, em, L)] = el
    phiE_lm = - 2 * E_lm / np.sqrt((ellfac+2)*(ellfac+1)*(ellfac)*(ellfac-1))
    phiE_lm[np.isnan(phiE_lm)] = 0
    kappa_lm = - ellfac * (ellfac + 1) * phiE_lm / 2.
    kappa_lm[np.isnan(kappa_lm)] = 0
    kappa_map_hp_rec = hp.alm2map(kappa_lm, nside=Nside, lmax=lmax, pol=False)

    return kappa_map_hp_rec


def gamma_to_kappa_hp(np.ndarray[double, ndim=1, mode="c"] gamma_real not None, \
    np.ndarray[double, ndim=1, mode="c"] gamma_imag not None, int L, int Nside, float sigma=-1):

    cdef int lmax=L-1

    maps_hp = [np.zeros(hp.nside2npix(Nside)),gamma_real, gamma_imag]
    [dummy, gamma_E_lm_rec, gamma_B_lm_rec] = hp.map2alm(maps_hp, lmax=lmax)

    kappa_E_lm, kappa_B_lm  = mw_mm.gamma_lm_to_kappa_lm_hp(gamma_E_lm_rec, gamma_B_lm_rec, L, sigma=sigma)

    kappa_map_E_hp = hp.alm2map(kappa_E_lm, nside=Nside, lmax=lmax, pol=False)
    kappa_map_B_hp = hp.alm2map(kappa_B_lm, nside=Nside, lmax=lmax, pol=False)

    return kappa_map_E_hp, kappa_map_B_hp

def reduced_shear_to_kappa_hp(np.ndarray[double, ndim=1, mode="c"] gamma_real not None, \
    np.ndarray[double, ndim=1, mode="c"] gamma_imag not None, int L, int Nside, float sigma=-1, \
    float tol_error=1E-10, bint Iterate=True, bint return_count=False):

    cdef np.ndarray[complex, ndim=1] kappa_E_lm, kappa_B_lm
    cdef np.ndarray[double, ndim=1] gamma_real_smooth, gamma_imag_smooth, gamma_real_dum, gamma_imag_dum
    cdef np.ndarray[double, ndim=1] kappa_E_1, kappa_E_2, kappa_B
    cdef int lmax=L-1, count, Npix=gamma_real.shape[0]
    cdef bint rel_error=True

    gamma_real_dum = np.zeros(Npix)
    gamma_imag_dum = np.zeros(Npix)

    if sigma>0:
        maps_hp = [np.zeros(hp.nside2npix(Nside)),gamma_real, gamma_imag]
        [dummy, gamma_real_smooth, gamma_imag_smooth] = hp.sphtfunc.smoothing(maps_hp, sigma=sigma, pol=True, lmax=lmax)
    else:
        gamma_real_smooth = gamma_real
        gamma_imag_smooth = gamma_imag

    maps_hp = [np.zeros(hp.nside2npix(Nside)),gamma_real_smooth, gamma_imag_smooth]
    [dummy, gamma_E_lm_rec, gamma_B_lm_rec] = hp.map2alm(maps_hp, lmax=lmax)

    kappa_E_lm, kappa_B_lm  = mw_mm.gamma_lm_to_kappa_lm_hp(gamma_E_lm_rec, gamma_B_lm_rec, L, sigma=-1)

    kappa_E_1 = hp.alm2map(kappa_E_lm, nside=Nside, lmax=lmax, pol=False)

    if Iterate:
        count=0
        while rel_error:
            count += 1
            if count>500:
                raise RuntimeError('Max iterations reached')

            for i in range(Npix):
                if not gamma_real[i] == hp.UNSEEN:
                    gamma_real_dum[i] = gamma_real_smooth[i]*(1-kappa_E_1[i].real)
                    gamma_imag_dum[i] = gamma_imag_smooth[i]*(1-kappa_E_1[i].real)

    # put transofrm in here
            maps_hp = [np.zeros(hp.nside2npix(Nside)),gamma_real_dum, gamma_imag_dum]
            [dummy, gamma_E_lm_rec, gamma_B_lm_rec] = hp.map2alm(maps_hp, lmax=lmax)

            kappa_E_lm, kappa_B_lm  = mw_mm.gamma_lm_to_kappa_lm_hp(gamma_E_lm_rec, gamma_B_lm_rec, L, sigma=-1)

            kappa_E_2 = hp.alm2map(kappa_E_lm, nside=Nside, lmax=lmax, pol=False)


            rel_error = False
            for i in range(Npix):
                if np.abs(kappa_E_2[i] -kappa_E_1[i]) > tol_error:
                        rel_error = True

            for i in range(Npix):
                kappa_E_1[i] = kappa_E_2[i]

# calculate kappa_b

    kappa_B = hp.alm2map(kappa_B_lm, nside=Nside, lmax=lmax, pol=False)

    if return_count:
        return kappa_E_1, kappa_B, count
    else:
        return kappa_E_1, kappa_B












    