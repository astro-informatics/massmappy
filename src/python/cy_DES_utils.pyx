import numpy as np
cimport numpy as np
import healpy as hp
import pyssht as ssht

#in this case makes a map of e1, correcting like Jarvis et al. Section 9.2
def make_shear_hp_map(np.ndarray[double, ndim=2, mode="c"] e1, np.ndarray[double, ndim=2, mode="c"] e2, \
	np.ndarray[double, ndim=2, mode="c"] c1, np.ndarray[double, ndim=2, mode="c"] c2, \
	np.ndarray[double, ndim=2, mode="c"] weight, \
	np.ndarray[double, ndim=2, mode="c"] mcorr, np.ndarray[long, ndim=2, mode="c"] pixnum, int Npix, int n_min):

	cdef int i, pix, dec_ring, RA_ring, Nside
	cdef np.ndarray[double, ndim=1] e1map, e2map, weightmap, mask
	cdef np.ndarray[long, ndim=1] Ngal

	e1map     = np.zeros(Npix)
	e2map     = np.zeros(Npix)
	weightmap = np.zeros(Npix)
	Ngal      = np.zeros(Npix, dtype=int)
	mask      = np.full((Npix,),np.nan)

	Nside = hp.npix2nside(Npix)

	for i in range(pixnum.size):
		pix = pixnum[i]
		e1map[pix]     += (weight[i,0]*(e1[i,0]-c1[i,0]))
		e2map[pix]     += (weight[i,0]*(e2[i,0]-c2[i,0]))#*(-1)
		weightmap[pix] += weight[i,0]*(1+mcorr[i,0])
		Ngal[pix]      += 1
		mask[pix]       = 1.0
#		print  pix, weightmap[pix], weight[i,0], (mcorr[i,0])

	for i in range(Npix):
		if Ngal[i]>n_min:
#			print i
			e1map[i] = e1map[i]/weightmap[i] 
			e2map[i] = e2map[i]/weightmap[i] 
		else:
			e1map[i] = hp.UNSEEN
			e2map[i] = hp.UNSEEN

	return e1map, e2map, mask

def make_shear_sine_map(np.ndarray[double, ndim=2, mode="c"] e1, np.ndarray[double, ndim=2, mode="c"] e2, \
	np.ndarray[double, ndim=2, mode="c"] c1, np.ndarray[double, ndim=2, mode="c"] c2, \
	np.ndarray[double, ndim=2, mode="c"] weight, \
	np.ndarray[double, ndim=2, mode="c"] mcorr, np.ndarray[double, ndim=2, mode="c"] RA, \
	np.ndarray[double, ndim=2, mode="c"] dec, double pix_size, int n_min, int Nx=215, int Ny=215, bint apply_rotation=True):

	cdef int i, j, pix_i, pix_j
	cdef float min_dec, pix_step_deg, psi_rot, rms_psi_rot=0.0
	cdef np.ndarray[double, ndim=2] e1map, e2map, weightmap, mask
	cdef np.ndarray[long, ndim=2] Ngal

	e1map     = np.zeros((Nx, Ny),)
	e2map     = np.zeros((Nx, Ny),)
	weightmap = np.zeros((Nx, Ny),)
	Ngal      = np.zeros((Nx, Ny), dtype=int)
	mask      = np.zeros((Nx, Ny),)

	pix_step_deg = pix_size/60.

	min_dec = dec.min()

	for i in range(dec.size):
		pix_j = int(70+np.cos(dec[i]*np.pi/180.)*(RA[i]-71)/pix_step_deg)
		pix_i = Nx-int((dec[i]-min_dec)/pix_step_deg)-1
#		print RA[i], dec[i], pix_i, pix_j

		if apply_rotation:
			psi_rot = -2*(RA[i]-71)*np.sin(dec[i]*np.pi/180)*np.pi/180
			e1map[pix_i,pix_j]     += weight[i,0]*((e1[i,0]-c1[i,0])*np.cos(psi_rot)-(e2[i,0]-c2[i,0])*np.sin(psi_rot))
			e2map[pix_i,pix_j]     += (weight[i,0]*(e2[i,0]-c2[i,0])*np.cos(psi_rot)+(e1[i,0]-c1[i,0])*np.sin(psi_rot))#*(-1)
		else:
			e1map[pix_i,pix_j]     += (weight[i,0]*(e1[i,0]-c1[i,0]))
			e2map[pix_i,pix_j]     += (weight[i,0]*(e2[i,0]-c2[i,0]))#*(-1)
		weightmap[pix_i,pix_j] += weight[i,0]*(1+mcorr[i,0])
		Ngal[pix_i,pix_j]      += 1


	for pix_i in range(Nx):
		for pix_j in range(Ny):
			if Ngal[pix_i,pix_j]>n_min:
#			print i
				e1map[pix_i,pix_j] = e1map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
				e2map[pix_i,pix_j] = e2map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
				mask[pix_i,pix_j]  = np.nan
			else:
				e1map[pix_i,pix_j] = np.nan
				e2map[pix_i,pix_j] = np.nan

	return e1map, e2map, mask


def make_kappa_sine_projection(np.ndarray[complex, ndim=2, mode="c"] k_mw, \
	double min_dec, int L, double pix_size=5.0, int Nx=215, \
	int Ny=215):

	cdef int i, j, pix_i, pix_j, n_theta, n_phi
	cdef float pix_step_deg, psi_rot, rms_psi_rot=0.0
	cdef np.ndarray[double, ndim=2] n_points, mask, theta, phi, dec, RA
	cdef np.ndarray[complex, ndim=2] k_sine

	k_sine   = np.zeros((Nx, Ny),dtype=complex)
	n_points = np.zeros((Nx, Ny),)
	mask     = np.zeros((Nx, Ny),)

	pix_step_deg = pix_size/60.

	n_theta, n_phi = ssht.sample_shape(L)
	theta, phi     = ssht.sample_positions(L, Grid=True)
	dec, RA        = ssht.theta_phi_to_ra_dec(theta,phi,Degrees=True)

	RA = (-1)*(RA-71.) + 71.

	for i in range(n_theta):
		for j in range(n_phi):
			if not np.isnan(k_mw[i,j]):
				pix_j = int(70+np.cos(dec[i,j]*np.pi/180.)*(RA[i,j]-71)/pix_step_deg)
				pix_i = Nx-int((dec[i,j]-min_dec)/pix_step_deg)-1
#				print RA[i,j], dec[i,j], min_dec, pix_i, pix_j

				k_sine[pix_i,pix_j]      = k_sine[pix_i,pix_j] + k_mw[i,j] 
				n_points[pix_i,pix_j]   += 1


	for pix_i in range(Nx):
		for pix_j in range(Ny):
			if n_points[pix_i,pix_j]>0:
				k_sine[pix_i,pix_j] = k_sine[pix_i,pix_j]/n_points[pix_i,pix_j]
				mask[pix_i,pix_j]   = np.nan
			else:
				k_sine[pix_i,pix_j] = np.nan

	return k_sine, mask



def make_shear_rot_sterio_map(np.ndarray[double, ndim=2, mode="c"] e1, np.ndarray[double, ndim=2, mode="c"] e2, \
	np.ndarray[double, ndim=2, mode="c"] c1, np.ndarray[double, ndim=2, mode="c"] c2, \
	np.ndarray[double, ndim=2, mode="c"] weight, \
	np.ndarray[double, ndim=2, mode="c"] mcorr, np.ndarray[double, ndim=2, mode="c"] RA, \
	np.ndarray[double, ndim=2, mode="c"] dec, double pix_size, int n_min, int Nx=250, int Ny=250):

	cdef int i, j, pix_i, pix_j
	cdef float pix_step_rad,   rho_pos,   x_pos,   y_pos, delta_x, delta_y, delta_z
	cdef float theta_pos, phi_pos, rho, delta_rho, delta_phi, delta_x_plane, delta_y_plane, psi_rot
	cdef np.ndarray[double, ndim=2] e1map, e2map, weightmap, mask, rot_matix
	cdef np.ndarray[double, ndim=2] theta,   phi,   xx,   yy,   zz,   xx_p,   yy_p,   zz_p
	cdef np.ndarray[long, ndim=2] Ngal

	e1map     = np.zeros((Nx, Ny),)
	e2map     = np.zeros((Nx, Ny),)
	weightmap = np.zeros((Nx, Ny),)
	Ngal      = np.zeros((Nx, Ny), dtype=int)
	mask      = np.zeros((Nx, Ny),)

	theta, phi =  ssht.ra_dec_to_theta_phi(RA, dec, Degrees=True)

	alpha = -21+180
	beta  = -37
	g     = 90

	rot = [np.radians(alpha),np.radians(beta),np.radians(g)]

	xx, yy, zz = ssht.s2_to_cart(theta, phi)
	xx_p, yy_p, zz_p = ssht.rot_cart_2d(xx, yy, zz, rot)
	theta, phi = ssht.cart_to_s2(xx_p, yy_p, zz_p)

	pix_step_rad = pix_size*np.pi/(60.*180.)

	rot_matix = ssht.make_rotation_matrix(rot)
	delta_x = rot_matix[0,2]
	delta_y = rot_matix[1,2]
	delta_z = rot_matix[2,2]

	for i in range(theta.size):
		theta_pos = theta[i]
		phi_pos   = phi[i]
		rho       = 2.0*np.tan((np.pi-theta_pos)/2.0)

		x_pos = rho*np.cos(phi_pos)
		y_pos = rho*np.sin(phi_pos)

		delta_rho  = np.sin(phi_pos)*delta_y+np.cos(phi_pos)*delta_x - np.tan(theta_pos)*delta_z
		delta_rho *= np.cos(theta_pos)
		delta_rho *= -1.0/(np.cos((np.pi-theta_pos)/2.0)*np.cos((np.pi-theta_pos)/2.0))

		delta_phi  = -np.tan(phi_pos)*delta_x + delta_y
		delta_phi *= np.cos(phi_pos)/np.sin(theta_pos)

		delta_x_plane = np.cos(phi_pos)*delta_rho - rho*np.sin(phi_pos)*delta_phi
		delta_y_plane = np.sin(phi_pos)*delta_rho + rho*np.cos(phi_pos)*delta_phi

		# calculate angle
		psi_rot = -2.0*np.arctan2(delta_x_plane, delta_y_plane)

		pix_i = int((x_pos)/pix_step_rad-0.5)+Nx/2
		pix_j = int((y_pos)/pix_step_rad-0.5)+Ny/2

#		e1map[pix_i,pix_j]     += weight[i,0]*psi_rot
#		psi_rot = -2*(RA[i]-71)*np.sin(dec[i]*np.pi/180)*np.pi/180
#		e2map[pix_i,pix_j]     += weight[i,0]*psi_rot
#		weightmap[pix_i,pix_j] += weight[i,0]
		e1map[pix_i,pix_j]     += weight[i,0]*((e1[i,0]-c1[i,0])*np.cos(psi_rot)-(e2[i,0]-c2[i,0])*np.sin(psi_rot))
		e2map[pix_i,pix_j]     += weight[i,0]*(e2[i,0]-c2[i,0])*np.cos(psi_rot)+(e1[i,0]-c1[i,0])*np.sin(psi_rot)
		weightmap[pix_i,pix_j] += weight[i,0]*(1+mcorr[i,0])
		Ngal[pix_i,pix_j]      += 1


	for pix_i in range(Nx):
		for pix_j in range(Ny):
			if Ngal[pix_i,pix_j]>n_min:
#			print i
				e1map[pix_i,pix_j] = e1map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
				e2map[pix_i,pix_j] = e2map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
				mask[pix_i,pix_j]  = np.nan
			else:
				e1map[pix_i,pix_j] = np.nan
				e2map[pix_i,pix_j] = np.nan

	return e1map, e2map, mask


def make_shear_mw_map(np.ndarray[double, ndim=2, mode="c"] e1, np.ndarray[double, ndim=2, mode="c"] e2, \
	np.ndarray[double, ndim=2, mode="c"] c1, np.ndarray[double, ndim=2, mode="c"] c2, \
	np.ndarray[double, ndim=2, mode="c"] weight, \
	np.ndarray[double, ndim=2, mode="c"] mcorr, np.ndarray[double, ndim=2, mode="c"] RA, \
	np.ndarray[double, ndim=2, mode="c"] dec, int n_min, int L=2160, str Method="MW"):

	cdef int i, j, pix_i, pix_j, dec_ring, RA_ring
	cdef float theta_i, phi_i
	cdef np.ndarray[double, ndim=2] e1map, e2map, weightmap
	cdef np.ndarray[double, ndim=2] theta, phi
	cdef np.ndarray[long, ndim=2] Ngal

	e1map     = np.zeros(ssht.sample_shape(L,Method=Method,),)
	e2map     = np.zeros(ssht.sample_shape(L,Method=Method,),)
	weightmap = np.zeros(ssht.sample_shape(L,Method=Method,),)
	Ngal      = np.zeros(ssht.sample_shape(L,Method=Method,), dtype=int)

	theta, phi =  ssht.ra_dec_to_theta_phi(RA, dec, Degrees=True)


	for i in range(theta.size):
		pix_i = ssht.theta_to_index(theta[i], L, Method=Method)
		pix_j = ssht.phi_to_index(phi[i], L, Method=Method)

		e1map[pix_i,pix_j]     += (weight[i,0]*(e1[i,0]-c1[i,0]))
		e2map[pix_i,pix_j]     += (weight[i,0]*(e2[i,0]-c2[i,0]))#*(-1)
		weightmap[pix_i,pix_j] += weight[i,0]*(1+mcorr[i,0])
		Ngal[pix_i,pix_j]      += 1
#		mask[pix_i,pix_j]       = 1.0


	for pix_i in range(e1map.shape[0]):
		for pix_j in range(e1map.shape[1]):
			if Ngal[pix_i,pix_j]>n_min:
#			print i
				e1map[pix_i,pix_j] = e1map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
				e2map[pix_i,pix_j] = e2map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
			else:
				e1map[pix_i,pix_j] = np.nan
				e2map[pix_i,pix_j] = np.nan


	return e1map, e2map




# def make_shear_rot_cyl_map(np.ndarray[double, ndim=2, mode="c"] e1, np.ndarray[double, ndim=2, mode="c"] e2, \
# 	np.ndarray[double, ndim=2, mode="c"] c1, np.ndarray[double, ndim=2, mode="c"] c2, \
# 	np.ndarray[double, ndim=2, mode="c"] weight, \
# 	np.ndarray[double, ndim=2, mode="c"] mcorr, np.ndarray[double, ndim=2, mode="c"] RA, \
# 	np.ndarray[double, ndim=2, mode="c"] dec, double pix_size, int n_min, int Nx=250, int Ny=200):

# 	cdef int i, j, pix_i, pix_j
# 	cdef float min_theta, min_phi, pix_step_rad
# 	cdef np.ndarray[double, ndim=2] e1map, e2map, weightmap, mask
# 	cdef np.ndarray[double, ndim=2] theta, phi, xx, yy, zz, xx_p, yy_p, zz_p
# 	cdef np.ndarray[long, ndim=2] Ngal

# 	e1map     = np.zeros((Nx, Ny),)
# 	e2map     = np.zeros((Nx, Ny),)
# 	weightmap = np.zeros((Nx, Ny),)
# 	Ngal      = np.zeros((Nx, Ny), dtype=int)
# 	mask      = np.full((Nx, Ny),np.nan)

# 	theta, phi =  ssht.ra_dec_to_theta_phi(RA, dec, Degrees=True)

# 	rot = [71*np.pi/180,np.pi/2,0.0]

# 	xx, yy, zz = ssht.s2_to_cart(theta, phi)
# 	xx_p, yy_p, zz_p = ssht.rot_cart_2d(xx, yy, zz, rot)
# 	theta, phi = ssht.cart_to_s2(xx_p, yy_p, zz_p)

# 	min_theta = theta.min()
# 	min_phi   = phi.min()
# 	pix_step_rad = pix_size*np.pi/(60.*180.)

# 	for i in range(theta.size):
# 		pix_i = int((phi[i]-min_phi)/pix_step_rad)
# 		pix_j = int((theta[i]-min_theta)/pix_step_rad)

# 		e1map[pix_i,pix_j]     += (weight[i,0]*(e1[i,0]-c1[i,0]))
# 		e2map[pix_i,pix_j]     += (weight[i,0]*(e2[i,0]-c2[i,0]))#*(-1)
# 		weightmap[pix_i,pix_j] += weight[i,0]*(1+mcorr[i,0])
# 		Ngal[pix_i,pix_j]      += 1
# 		mask[pix_i,pix_j]       = 1.0


# 	for pix_i in range(Nx):
# 		for pix_j in range(Ny):
# 			if Ngal[pix_i,pix_j]>n_min:
# #			print i
# 				e1map[pix_i,pix_j] = e1map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
# 				e2map[pix_i,pix_j] = e2map[pix_i,pix_j]/weightmap[pix_i,pix_j] 
# 			else:
# 				e1map[pix_i,pix_j] = np.nan
# 				e2map[pix_i,pix_j] = np.nan

# 	return e1map, e2map, mask