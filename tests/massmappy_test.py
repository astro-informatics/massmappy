import numpy as np
import cy_mass_mapping as mm
import cy_healpy_mass_mapping as hp_mm
import pyssht as ssht


def test_lm_conversion_functions():
	''' To test the lm conversion functions
	'''

	L = 3

	flm_mw = np.array([0., -1+1j, 1+0j, 1+1j, 2-2j, -2+1j, 2+0j, 2+1j, 2+2j])
	flm_hp = np.array([0., 1+0j, 2+0j, 1+1j, 2+1j, 2+2j])
	
	flm_hp_test = mm.lm2lm_hp(flm_mw, 3)
	flm_mw_test = mm.lm_hp2lm(flm_hp, 3)

	np.testing.assert_almost_equal(flm_hp_test, flm_hp, decimal=7, err_msg='mm.lm2lm_hp failed to convert properly', verbose=True)
	np.testing.assert_almost_equal(flm_mw_test, flm_mw, decimal=7, err_msg='mm.lm2lm_mw failed to convert properly', verbose=True)


def test_hp_index2lm():
	''' To test the MW index functions
	'''

	L = 3


	np.testing.assert_equal(mm.healpy_ind2lm(0,L)[0],0,err_msg='hp_mm.healpy_ind2lm failed to get the correct ell')
	np.testing.assert_equal(mm.healpy_ind2lm(0,L)[1],0,err_msg='hp_mm.healpy_ind2lm failed to get the correct em')

	np.testing.assert_equal(mm.healpy_ind2lm(1,L)[0],1,err_msg='hp_mm.healpy_ind2lm failed to get the correct ell')
	np.testing.assert_equal(mm.healpy_ind2lm(1,L)[1],0,err_msg='hp_mm.healpy_ind2lm failed to get the correct em')

	np.testing.assert_equal(mm.healpy_ind2lm(2,L)[0],2,err_msg='hp_mm.healpy_ind2lm failed to get the correct ell')
	np.testing.assert_equal(mm.healpy_ind2lm(2,L)[1],0,err_msg='hp_mm.healpy_ind2lm failed to get the correct em')

	np.testing.assert_equal(mm.healpy_ind2lm(3,L)[0],1,err_msg='hp_mm.healpy_ind2lm failed to get the correct ell')
	np.testing.assert_equal(mm.healpy_ind2lm(3,L)[1],1,err_msg='hp_mm.healpy_ind2lm failed to get the correct em')

	np.testing.assert_equal(mm.healpy_ind2lm(4,L)[0],2,err_msg='hp_mm.healpy_ind2lm failed to get the correct ell')
	np.testing.assert_equal(mm.healpy_ind2lm(4,L)[1],1,err_msg='hp_mm.healpy_ind2lm failed to get the correct em')

	np.testing.assert_equal(mm.healpy_ind2lm(5,L)[0],2,err_msg='hp_mm.healpy_ind2lm failed to get the correct ell')
	np.testing.assert_equal(mm.healpy_ind2lm(5,L)[1],2,err_msg='hp_mm.healpy_ind2lm failed to get the correct em')

def test_hp_lm2index():

	L=3

	np.testing.assert_equal(mm.healpy_lm2ind(0,0,L),0,err_msg='hp_mm.healpy_lm2ind failed to get the correct index')

	np.testing.assert_equal(mm.healpy_lm2ind(1,0,L),1,err_msg='hp_mm.healpy_lm2ind failed to get the correct index')

	np.testing.assert_equal(mm.healpy_lm2ind(2,0,L),2,err_msg='hp_mm.healpy_lm2ind failed to get the correct index')

	np.testing.assert_equal(mm.healpy_lm2ind(1,1,L),3,err_msg='hp_mm.healpy_lm2ind failed to get the correct index')

	np.testing.assert_equal(mm.healpy_lm2ind(2,1,L),4,err_msg='hp_mm.healpy_lm2ind failed to get the correct index')

	np.testing.assert_equal(mm.healpy_lm2ind(2,2,L),5,err_msg='hp_mm.healpy_lm2ind failed to get the correct index')

def test_generate_kappa_lm_hp():

	L = 3
	Cl = np.array([1E2,1E4,1E6])
	seed = 1
	mm.generate_kappa_lm_hp(Cl, L, seed=seed)

	kappa_prep = np.array([0.0+0j, 0.0+0j,1624.34536366+0j, 0.0+0j,\
        -611.75641365-528.17175226j,-1072.96862216+865.40762932j])

	kappa_gen = mm.generate_kappa_lm_hp(Cl, L, seed=seed)

	np.testing.assert_almost_equal(kappa_gen, kappa_prep, decimal=7, err_msg='mm.generate_kappa_lm_hp failed to generate standard output', verbose=True)


def test_generate_kappa_lm_mw():

	L = 3
	Cl = np.array([1E2,1E4,1E6])
	seed = 1

	kappa_gen = mm.generate_kappa_lm_mw(Cl, L, seed=seed)

	kappa_prep = np.array([0.0+0j, 0.0+0j, 0.0+0j, 0.0+0j,\
       -1072.96862216-865.40762932j,   611.75641365-528.17175226j,\
        1624.34536366  +0j,  -611.75641365-528.17175226j,\
       -1072.96862216+865.40762932j])

	np.testing.assert_almost_equal(kappa_gen, kappa_prep, decimal=7, err_msg='mm.generate_kappa_lm_hp failed to generate standard output', verbose=True)

def test_lm_transform_mw():

	L = 20
	Cl = np.ones(L)
	seed = 1

	kappa_lm = mm.generate_kappa_lm_mw(Cl, L, seed=seed)

	gamma_lm = mm.kappa_lm_to_gamma_lm_mw(kappa_lm, L)

	kappa_lm_rec = mm.gamma_lm_to_kappa_lm_mw(gamma_lm, L)

	np.testing.assert_almost_equal(kappa_lm_rec, kappa_lm, decimal=7, err_msg="lm gamma and kappa MW transoforms are not invereses")

def test_lm_transform_hp():

	L = 20
	Cl = np.ones(L)
	seed = 1

	kappa_lm = mm.generate_kappa_lm_hp(Cl, L, seed=seed)

	gamma_lm_E, gamma_lm_B = mm.kappa_lm_to_gamma_lm_hp(kappa_lm, L)

	kappa_lm_E_rec, kappa_lm_B_rec = mm.gamma_lm_to_kappa_lm_hp(gamma_lm_E, gamma_lm_B, L)

	np.testing.assert_almost_equal(kappa_lm_E_rec, kappa_lm, decimal=7, err_msg="lm gamma and kappa HEALPix transoforms are not invereses (E)")
	np.testing.assert_almost_equal(kappa_lm_B_rec+1.0, np.ones(kappa_lm.size), decimal=7, err_msg="lm gamma and kappa HEALPix transoforms are not invereses (B)")

def test_gamma_transform_mw():

	L = 20
	Cl = np.ones(L)
	seed = 1
	Method = "MW"

	kappa_lm = mm.generate_kappa_lm_mw(Cl, L, seed=seed)

	k_mw = ssht.inverse(kappa_lm, L, Reality=True, Method=Method)

	gamma_lm = mm.kappa_lm_to_gamma_lm_mw(kappa_lm, L)

	gamma = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)

	k_rec_mw = mm.gamma_to_kappa_mw(gamma, L, Method=Method)

	np.testing.assert_almost_equal(k_rec_mw, k_mw, decimal=7, err_msg="gamma and kappa transoforms are not invereses")

def test_reduced_shear_transform_mw():

	L = 20
	Cl = np.ones(L)*1E-4
	seed = 1
	Method = "MW"

	kappa_lm = mm.generate_kappa_lm_mw(Cl, L, seed=seed)

	k_mw = ssht.inverse(kappa_lm, L, Reality=True, Method=Method)

	gamma_lm = mm.kappa_lm_to_gamma_lm_mw(kappa_lm, L)

	gamma = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)

	shear = gamma/(1.0-k_mw)

	k_rec_mw = mm.reduced_shear_to_kappa_mw(shear, L, Iterate=True, Method=Method, tol_error=1E-10)

	np.testing.assert_almost_equal(k_rec_mw, k_mw, decimal=7, err_msg="reduced shear and kappa transoforms are not invereses")


if __name__ == '__main__':
	test_lm_conversion_functions()     #1
	test_hp_index2lm()                 #2
	test_hp_lm2index()                 #3
	test_generate_kappa_lm_hp()        #4
	test_generate_kappa_lm_mw()        #5
	test_lm_transform_mw()             #6
	test_lm_transform_hp()             #7
	test_gamma_transform_mw()          #8
	test_reduced_shear_transform_mw()  #9