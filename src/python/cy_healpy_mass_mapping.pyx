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