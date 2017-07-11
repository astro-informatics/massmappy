import numpy as np
cimport numpy as np
import healpy as hp
import pyssht as ssht

from astropy.io import fits
import scipy.io as sio
from matplotlib import pyplot as plt


def make_shear_map_hp(filenames, int N_side, int N_mingal=10):

	cdef np.ndarray[double, ndim=1] e1map, e2map, weight, RA, DEC, e1, e2
	cdef np.ndarray[np.int64_t, ndim=1] pixnum

	cdef int i_fits, i_map, N_pix=12*N_side*N_side, N_fits, i_file

	e1map  = np.zeros(N_pix)
	e2map  = np.zeros(N_pix)
	weight = np.zeros(N_pix)

	if isinstance(filenames, basestring):
		filenames = [filenames]

	for i_file, filename in enumerate(filenames):
		cat = fits.open(filename)

		RA  = cat[1].data["ra_arcmin"]*np.pi/(180)
		DEC = cat[1].data["dec_arcmin"]*np.pi/(180)
		e1  = cat[1].data["shear1"]*1.0  # don't ask! its a bad work around due to some issue with astropy (according to some guy on stack exachange)
		e2  = cat[1].data["shear2"]*1.0

		N_fits = RA.size
		print "In file " + str(i_file) + " we have: ", N_fits, " galaxies"

		DEC = np.pi/2 - DEC
		pixnum = hp.ang2pix(N_side, DEC, RA)

		for i_fits in range(N_fits):
			e1map[pixnum[i_fits]]  += e1[i_fits]
			e2map[pixnum[i_fits]]  += e2[i_fits]
			weight[pixnum[i_fits]] += 1

	e1map[weight>N_mingal-1] = e1map[weight>N_mingal-1]/weight[weight>N_mingal-1]
	e2map[weight>N_mingal-1] = e2map[weight>N_mingal-1]/weight[weight>N_mingal-1]

	e1map[weight<N_mingal]  = hp.UNSEEN
	e2map[weight<N_mingal]  = hp.UNSEEN

	# hp.mollview(e1map)
	# hp.mollview(e2map)
	# hp.mollview(weight)
	# plt.show()

	return e1map, e2map, weight


def make_shear_map_mw(filenames, int L, int N_mingal=10):

	cdef np.ndarray[double, ndim=1] RA, DEC, e1, e2, x, y, z, x_p, y_p, z_p
	cdef np.ndarray[double, ndim=2] e1map, e2map, weight
	cdef np.ndarray[np.int64_t, ndim=1] pixnum

	cdef int i_fits, i_map, N_theta, N_phi, N_fits, i_file, i_theta, i_phi

	N_theta, N_phi = ssht.sample_shape(L)

	e1map  = np.zeros((N_theta,N_phi))
	e2map  = np.zeros((N_theta,N_phi))
	weight = np.zeros((N_theta,N_phi))

	if isinstance(filenames, basestring):
		filenames = [filenames]

	for i_file, filename in enumerate(filenames):
		cat = fits.open(filename)

		RA  = cat[1].data["ra_arcmin"]*np.pi/(180)
		DEC = cat[1].data["dec_arcmin"]*np.pi/(180)
		e1  = cat[1].data["shear1"]*1.0  # don't ask! its a bad work around due to some issue with astropy (according to some guy on stack exachange)
		e2  = cat[1].data["shear2"]*1.0

		N_fits = RA.size
		print "In file " + str(i_file) + " we have: ", N_fits, " galaxies"

		DEC = np.pi/2 - DEC

		# x, y, z = ssht.s2_to_cart(DEC, RA)
		# x_p, y_p, z_p = ssht.rot_cart_1d(x, y, z, [0.0, np.pi/4, np.pi/2])
		# DEC, RA = ssht.cart_to_s2(x_p, y_p, z_p)

		for i_fits in range(N_fits):
			i_theta =  ssht.theta_to_index(DEC[i_fits],L)
			i_phi   =  ssht.phi_to_index(RA[i_fits],L)
			e1map[i_theta, i_phi]  += e1[i_fits]
			e2map[i_theta, i_phi]  += e2[i_fits]
			weight[i_theta, i_phi] += 1

	e1map[weight>N_mingal-1] = e1map[weight>N_mingal-1]/weight[weight>N_mingal-1]
	e2map[weight>N_mingal-1] = e2map[weight>N_mingal-1]/weight[weight>N_mingal-1]

	e1map[weight<N_mingal]  = np.nan
	e2map[weight<N_mingal]  = np.nan

	plt.imshow(weight)
	plt.show()

	return e1map, e2map, weight
