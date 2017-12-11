#!/usr/bin/python
'''
Read the DES data into Healpix
'''

#importing packages
import scipy.io as sio
import numpy as np
import scipy as sp
import scipy.ndimage as im
import random as rn
import healpy as hp
import matplotlib.mlab as mlab
import pyssht as ssht
import massmappy.cy_mass_mapping as mm
import massmappy.cy_healpy_mass_mapping as hp_mm
import massmappy.cy_DES_utils as DES
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import cm


def sterio_grid_lines(Projection="SP", RA_lines_points=[50, 55, 60, 65, 70, 75, 80, 85, 90, 95], dec_lines_points=[-45, -50, -55, -60], \
	rot=[np.radians(-21.+180.),np.radians(-37.),np.radians(90.)], resolution=500, zoom_region=np.pi/16*0.9):

	N = 200
	rot_inv = rot

	RA_lines_1  = np.zeros((N,len(RA_lines_points)), dtype=float)
	RA_lines_2  = np.zeros((N,len(RA_lines_points)), dtype=float)
	dec_lines_1 = np.zeros((N,len(dec_lines_points)), dtype=float)
	dec_lines_2 = np.zeros((N,len(dec_lines_points)), dtype=float)

	for i, point in enumerate(RA_lines_points):
		RA_lines_1[:,i] = (-1)*(point-71.) + 71.
		RA_lines_2[:,i] = np.linspace(-40,-70,N)

	for i, point in enumerate(dec_lines_points):
		dec_lines_1[:,i] = (-1)*(np.linspace(50,120,N)-71.) + 71.
		dec_lines_2[:,i] = point

	dec_lines_1 = dec_lines_1*np.pi/180
	dec_lines_2 = dec_lines_2*np.pi/180
	RA_lines_1  = RA_lines_1*np.pi/180
	RA_lines_2  = RA_lines_2*np.pi/180

	RA_lines_2[:,:]  = np.pi/2 - RA_lines_2[:,:]
	dec_lines_2[:,:] = np.pi/2 - dec_lines_2[:,:]

	xx, yy, zz = ssht.s2_to_cart(RA_lines_2, RA_lines_1)
	xx_p, yy_p, zz_p = ssht.rot_cart(xx, yy, zz, rot_inv)
	RA_lines_2, RA_lines_1 = ssht.cart_to_s2(xx_p, yy_p, zz_p)

	xx, yy, zz = ssht.s2_to_cart(dec_lines_2, dec_lines_1)
	xx_p, yy_p, zz_p = ssht.rot_cart_2d(xx, yy, zz, rot_inv)
	dec_lines_2, dec_lines_1 = ssht.cart_to_s2(xx_p, yy_p, zz_p)

	half_box_len = 2.0*np.tan(zoom_region/2)

	rho = 2.0*np.tan((np.pi-RA_lines_2)/2.0)
	phi = RA_lines_1.copy()


	RA_lines_1 = rho*np.cos(phi)
	RA_lines_2 = rho*np.sin(phi)

	RA_lines_1[:,:] = resolution*(RA_lines_1[:,:]+half_box_len)/(2.0*half_box_len)
	RA_lines_2[:,:] = resolution*(RA_lines_2[:,:]+half_box_len)/(2.0*half_box_len)

	rho = 2.0*np.tan((np.pi-dec_lines_2[:,:])/2.0)
	phi = dec_lines_1[:,:]

	dec_lines_1 = rho*np.cos(phi)
	dec_lines_2 = rho*np.sin(phi)

	dec_lines_1[:,:] = resolution*(dec_lines_1[:,:]+half_box_len)/(2.0*half_box_len)
	dec_lines_2[:,:] = resolution*(dec_lines_2[:,:]+half_box_len)/(2.0*half_box_len)

	# for i in range(len(RA_lines_points)):
	# 	for j in range (N):
	# 		if  np.abs(RA_lines_1[j,i]-resolution) < 5:
	# 			print RA_lines_1[j,i], RA_lines_2[j,i]
	# 	print ""

	# for i in range(len(dec_lines_points)):
	# 	for j in range (N):
	# 		if  np.abs(dec_lines_2[j,i]) < 5:
	# 			print dec_lines_2[j,i], dec_lines_1[j,i]
	# 	print ""


	return RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2

def sine_grid_lines(dec, RA_lines_points=[50, 55, 60, 65, 70, 75, 80, 85, 90, 95], dec_lines_points=[-45, -50, -55, -60], \
	resolution=215, pix_step_arcmin=5.0):

	N = 200

	RA_lines_1  = np.zeros((N,len(RA_lines_points)), dtype=float)
	RA_lines_2  = np.zeros((N,len(RA_lines_points)), dtype=float)
	dec_lines_1 = np.zeros((N,len(dec_lines_points)), dtype=float)
	dec_lines_2 = np.zeros((N,len(dec_lines_points)), dtype=float)

	for i, point in enumerate(RA_lines_points):
		RA_lines_1[:,i] = point
		RA_lines_2[:,i] = np.linspace(-40,-70,N)

	for i, point in enumerate(dec_lines_points):
		dec_lines_1[:,i] = np.linspace(50,120,N)
		dec_lines_2[:,i] = point



	min_dec = dec.min()
	pix_step_deg = pix_step_arcmin/60.

	RA_lines_1 = 70+np.cos(RA_lines_2*np.pi/180.)*(RA_lines_1-71)/pix_step_deg
	RA_lines_2 = resolution-(RA_lines_2-min_dec)/pix_step_deg-1


	dec_lines_1 = 70+np.cos(dec_lines_2*np.pi/180.)*(dec_lines_1-71)/pix_step_deg
	dec_lines_2 = resolution-(dec_lines_2-min_dec)/pix_step_deg-1

	# for i in range(len(RA_lines_points)):
	# 	for j in range (N):
	# 		if  np.abs(RA_lines_2[j,i]-resolution) < 5:
	# 			print RA_lines_2[j,i], RA_lines_1[j,i]
	# 	print ""

	# for i in range(len(dec_lines_points)):
	# 	for j in range (N):
	# 		if  np.abs(dec_lines_1[j,i]) < 5:
	# 			print dec_lines_1[j,i], dec_lines_2[j,i]
	# 	print ""


	return RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2

mat_contents = sio.loadmat('data/DES.mat')

Nside  = 512
L_hp = 4*Nside
L_mw = 2160 # 2160
Iterate = True

save_figs = True
show_figs = True


Npix   = hp.nside2npix(Nside)
ra     = np.ascontiguousarray(mat_contents['RA'])
ra_flip = (-1)*(ra-71.) + 71.
dec    = np.ascontiguousarray(mat_contents['dec'])
theta, phi = ssht.ra_dec_to_theta_phi(ra_flip, dec, Degrees=True)
pixnum = hp.ang2pix(Nside, theta, phi)

e1=np.ascontiguousarray(mat_contents['e1'])
e2=np.ascontiguousarray(mat_contents['e2'])
e2_flip = e2.copy()*(-1)

mcorr  = np.ascontiguousarray(mat_contents['mcorr'])
c1     = np.ascontiguousarray(mat_contents['c1'])
c2     = np.ascontiguousarray(mat_contents['c2'])
c2_flip = c2.copy()*(-1)
weight     = np.ascontiguousarray(mat_contents['weight'])

alpha = -21+180
beta  = -37
g     = 90

zoom_region  = 2.*np.arctan(0.25*250.*5.0*np.pi/(180.*60.))
sigma_scale = 1.0

e1map_hp, e2map_hp, mask_hp = DES.make_shear_hp_map(e1, e2, c1, c2, weight, mcorr, pixnum, Npix, 2)

e1map_hp[e1map_hp!=hp.UNSEEN] = e1map_hp[e1map_hp!=hp.UNSEEN]*(-1)
e2map_hp[e2map_hp!=hp.UNSEEN] = e2map_hp[e2map_hp!=hp.UNSEEN]*(-1)

mask = np.ones(Npix)
mask[e2map_hp==hp.UNSEEN] = hp.UNSEEN

# print Npix*float(ra.size)/(float(mask[mask!=hp.UNSEEN].size)*41253*60*60)

# hp.gnomview(e1map_hp, title="e1 map",rot=[70,-53,0.0], reso=6.0)
# hp.gnomview(e2map_hp, title="e2 map",rot=[70,-53,0.0], reso=6.0)


sigma = np.sqrt(2)*2.*20.0*np.pi/(60.*180.*2.355)

kappa_E_map_hp_rec, kappa_B_map_hp_rec = hp_mm.reduced_shear_to_kappa_hp(e1map_hp, e2map_hp, L_hp, Nside, sigma=sigma*sigma_scale, Iterate=Iterate)
# kappa_E_map_hp_rec, kappa_B_map_hp_rec = hp_mm.gamma_to_kappa_hp(e1map_hp, e2map_hp, L_hp, Nside, sigma=sigma*sigma_scale)


# hp.gnomview(kappa_E_map_hp_rec*mask_hp, title="K E map",rot=[70,-53,0.0], reso=6.0, min=-0.015, max=0.015, cmap="cubehelix")
# if save_figs:
# 	plt.savefig("fig/DES_hp_E.pdf")
# hp.gnomview(kappa_B_map_hp_rec*mask_hp, title="K B map",rot=[70,-53,0.0], reso=6.0, min=-0.015, max=0.015, cmap="cubehelix")
# if save_figs:
# 	plt.savefig("fig/DES_hp_B.pdf")

alm_mask_hp = hp.map2alm(mask, lmax=L_hp-1)
alm_E_hp    = hp.map2alm(kappa_E_map_hp_rec, lmax=L_hp-1)
alm_B_hp    = hp.map2alm(kappa_B_map_hp_rec, lmax=L_hp-1)

alm_mask_hp_mw = mm.lm_hp2lm(alm_mask_hp, L_hp)
alm_E_hp_mw    = mm.lm_hp2lm(alm_E_hp, L_hp)
alm_B_hp_mw    = mm.lm_hp2lm(alm_B_hp, L_hp)

mask_mw               = ssht.inverse(alm_mask_hp_mw, L_hp, Reality=True)
kappa_E_map_hp_rec_mw = ssht.inverse(alm_E_hp_mw, L_hp, Reality=True)
kappa_B_map_hp_rec_mw = ssht.inverse(alm_B_hp_mw, L_hp, Reality=True)

mask_mw[mask_mw<0.5] = np.nan

k_hp_mw = kappa_E_map_hp_rec_mw + 1j*kappa_B_map_hp_rec_mw

k_mw_north_real, mask_north_real, k_mw_south_real, mask_south_real, \
k_mw_north_imag, mask_north_imag, k_mw_south_imag, mask_south_imag \
		= ssht.polar_projection(k_hp_mw*mask_mw, L_hp, resolution=250, Method="MW", zoom_region=zoom_region,\
			rot=[np.radians(alpha),np.radians(beta),np.radians(g)], Projection="SP")

RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2 = sterio_grid_lines(resolution=250,zoom_region=zoom_region)

fig, ax = plt.subplots()
imgplot = ax.imshow(k_mw_south_real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
plt.colorbar(imgplot)
ax.imshow(mask_south_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(10):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(4):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
#plt.axis('off')
if save_figs:
	plt.savefig("fig/DES_hp_E_sterio.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(k_mw_south_imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
plt.colorbar(imgplot)
ax.imshow(mask_south_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(4):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
#plt.axis('off')
if save_figs:
	plt.savefig("fig/DES_hp_B_sterio.pdf")

e1map_hp, e2map_hp, mask_hp = DES.make_shear_hp_map(e1, e2, c1, c2, weight, mcorr, pixnum, Npix, 2)

e1map_hp[e1map_hp!=hp.UNSEEN] = e1map_hp[e1map_hp!=hp.UNSEEN]*(-1)
e2map_hp[e2map_hp!=hp.UNSEEN] = e2map_hp[e2map_hp!=hp.UNSEEN]*(-1)

mask = np.ones(Npix)
mask[e2map_hp==hp.UNSEEN] = hp.UNSEEN

print Npix*float(ra.size)/(float(mask[mask!=hp.UNSEEN].size)*41253*60*60)

# hp.gnomview(e1map_hp, title="e1 map",rot=[70,-53,0.0], reso=6.0)
# hp.gnomview(e2map_hp, title="e2 map",rot=[70,-53,0.0], reso=6.0)


sigma = np.sqrt(2)*2.*20.0*np.pi/(60.*180.*2.355)

kappa_E_map_hp_rec, kappa_B_map_hp_rec = hp_mm.reduced_shear_to_kappa_hp(e1map_hp, e2map_hp, L_hp, Nside, sigma=sigma*sigma_scale, Iterate=Iterate)
# kappa_E_map_hp_rec, kappa_B_map_hp_rec = hp_mm.gamma_to_kappa_hp(e1map_hp, e2map_hp, L_hp, Nside, sigma=sigma*sigma_scale)


# hp.gnomview(kappa_E_map_hp_rec*mask_hp, title="K E map",rot=[70,-53,0.0], reso=6.0, min=-0.015, max=0.015, cmap="cubehelix")
# if save_figs:
# 	plt.savefig("fig/DES_hp_E.pdf")
# hp.gnomview(kappa_B_map_hp_rec*mask_hp, title="K B map",rot=[70,-53,0.0], reso=6.0, min=-0.015, max=0.015, cmap="cubehelix")
# if save_figs:
# 	plt.savefig("fig/DES_hp_B.pdf")

alm_mask_hp = hp.map2alm(mask, lmax=L_hp-1)
alm_E_hp    = hp.map2alm(kappa_E_map_hp_rec, lmax=L_hp-1)
alm_B_hp    = hp.map2alm(kappa_B_map_hp_rec, lmax=L_hp-1)

alm_mask_hp_mw = mm.lm_hp2lm(alm_mask_hp, L_hp)
alm_E_hp_mw    = mm.lm_hp2lm(alm_E_hp, L_hp)
alm_B_hp_mw    = mm.lm_hp2lm(alm_B_hp, L_hp)

mask_mw               = ssht.inverse(alm_mask_hp_mw, L_hp, Reality=True)
kappa_E_map_hp_rec_mw = ssht.inverse(alm_E_hp_mw, L_hp, Reality=True)
kappa_B_map_hp_rec_mw = ssht.inverse(alm_B_hp_mw, L_hp, Reality=True)

mask_mw[mask_mw<0.5] = np.nan

k_hp_mw = kappa_E_map_hp_rec_mw + 1j*kappa_B_map_hp_rec_mw

k_mw_north_real, mask_north_real, k_mw_south_real, mask_south_real, \
k_mw_north_imag, mask_north_imag, k_mw_south_imag, mask_south_imag \
		= ssht.polar_projection(k_hp_mw*mask_mw, L_hp, resolution=250, Method="MW", zoom_region=zoom_region,\
			rot=[np.radians(alpha),np.radians(beta),np.radians(g)], Projection="SP")

RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2 = sterio_grid_lines(resolution=250,zoom_region=zoom_region)

fig, ax = plt.subplots()
imgplot = ax.imshow(k_mw_south_real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
plt.colorbar(imgplot)
ax.imshow(mask_south_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(10):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(4):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
#plt.axis('off')
if save_figs:
	plt.savefig("fig/DES_hp_E_sterio.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(k_mw_south_imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
plt.colorbar(imgplot)
ax.imshow(mask_south_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(4):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
#plt.axis('off')
if save_figs:
	plt.savefig("fig/DES_hp_B_sterio.pdf")

e1map_mw, e2map_mw= DES.make_shear_mw_map(e1, e2, c1, c2, weight, mcorr, ra_flip, dec, 2, L=L_mw, Method="MW")


gamma_mw =  e1map_mw + 1j*e2map_mw

k_mw = mm.reduced_shear_to_kappa_mw(gamma_mw, L_mw, Method="MW", sigma=2.*20.0*np.pi/(60.*180.*2.355)*sigma_scale, Iterate=Iterate)

k_mw_north_real, mask_north_real, k_mw_south_real, mask_south_real, \
k_mw_north_imag, mask_north_imag, k_mw_south_imag, mask_south_imag \
		= ssht.polar_projection(k_mw, L_mw, resolution=250, Method="MW", zoom_region=zoom_region,\
			rot=[np.radians(alpha),np.radians(beta),np.radians(g)], Projection="SP")

RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2 = sterio_grid_lines(resolution=250,zoom_region=zoom_region)

fig, ax = plt.subplots()
imgplot = ax.imshow(k_mw_south_real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
plt.colorbar(imgplot)
ax.imshow(mask_south_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(10):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(4):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
#plt.axis('off')
if save_figs:
	plt.savefig("fig/DES_MW_E.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(k_mw_south_imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
plt.colorbar(imgplot)
ax.imshow(mask_south_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(4):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
#plt.axis('off')
if save_figs:
	plt.savefig("fig/DES_MW_B.pdf")


# f_plot_real, masked_array_real, f_plot_imag, masked_array_imag = ssht.mollweide_projection(gamma_mw, L, resolution=500, Method="MW")
# plt.figure()
# imgplot = plt.imshow(f_plot_real,interpolation='nearest')
# plt.colorbar(imgplot,fraction=0.025, pad=0.04)
# plt.imshow(masked_array_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.gca().set_aspect("equal")
# plt.title("Sphere kappa")
# plt.axis('off')

# f_plot_real, masked_array_real, f_plot_imag, masked_array_imag = ssht.mollweide_projection(k_mw, L, resolution=500, Method="MW")
# plt.figure()
# imgplot = plt.imshow(f_plot_real,interpolation='nearest')
# plt.colorbar(imgplot,fraction=0.025, pad=0.04)
# plt.imshow(masked_array_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.gca().set_aspect("equal")
# plt.title("Sphere kappa")
# plt.axis('off')



e1map_sine, e2map_sine, mask_sine = DES.make_shear_sine_map(e1, e2, c1, c2, weight, mcorr, ra, dec, 5.0, 2, \
	Nx=215, Ny=215, apply_rotation=True)


# plt.figure()
# imgplot = plt.imshow(e1map_sine,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.title("e1_sine")
# plt.figure()
# imgplot = plt.imshow(e2map_sine,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.title("e2_sine")


gamma_plane = - e1map_sine + 1j*e2map_sine

kappa_sine = mm.reduced_shear_to_kappa_plane(gamma_plane, 1.0, 1.0, sigma=2*8./2.355*sigma_scale, Iterate=Iterate)

RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2 = sine_grid_lines(dec)

fig, ax = plt.subplots()
imgplot = ax.imshow(kappa_sine.real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
plt.colorbar(imgplot)
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
if save_figs:
	plt.savefig("fig/DES_sine_E.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(-kappa_sine.imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sine_B.pdf")

k_mw_sine, mask_mw_sine = DES.make_kappa_sine_projection(k_mw, dec.min(), L_mw)

fig, ax = plt.subplots()
#imgplot = ax.imshow(kappa_sine.real-k_mw_sine.real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
imgplot = ax.imshow(kappa_sine.real-k_mw_sine.real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.imshow(mask_mw_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
plt.colorbar(imgplot)
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
if save_figs:
	plt.savefig("fig/DES_sine_mw_diff_E.pdf")

fig, ax = plt.subplots()
#imgplot = ax.imshow(-kappa_sine.imag-k_mw_sine.imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
imgplot = ax.imshow(-kappa_sine.imag-k_mw_sine.imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.imshow(mask_mw_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sine_mw_diff_B.pdf")

e1map_sine_false, e2map_sine_false, mask_sine = DES.make_shear_sine_map(e1, e2, c1, c2, weight, mcorr, ra, dec, 5.0, 2, \
	Nx=215, Ny=215, apply_rotation=False)


# plt.figure()
# imgplot = plt.imshow(e1map_sine,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.title("e1_sine")
# plt.figure()
# imgplot = plt.imshow(e2map_sine,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.title("e2_sine")


gamma_plane = - e1map_sine_false + 1j*e2map_sine_false

kappa_sine_false = mm.reduced_shear_to_kappa_plane(gamma_plane, 1.0, 1.0, sigma=2*8./2.355*sigma_scale, Iterate=Iterate)

RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2 = sine_grid_lines(dec)

fig, ax = plt.subplots()
imgplot = ax.imshow(kappa_sine_false.real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
plt.colorbar(imgplot)
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
if save_figs:
	plt.savefig("fig/DES_sine_false_E.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(-kappa_sine_false.imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sine_false_B.pdf")

fig, ax = plt.subplots()
#imgplot = ax.imshow(kappa_sine.real-k_mw_sine.real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
imgplot = ax.imshow(kappa_sine.real-kappa_sine_false.real,interpolation='nearest',vmin=-0.005,vmax=0.005, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
plt.colorbar(imgplot)
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
if save_figs:
	plt.savefig("fig/DES_sine_diff_false_E.pdf")

fig, ax = plt.subplots()
#imgplot = ax.imshow(-kappa_sine.imag-k_mw_sine.imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
imgplot = ax.imshow(-kappa_sine.imag+kappa_sine_false.imag,interpolation='nearest',vmin=-0.005,vmax=0.005, cmap="cubehelix")
ax.imshow(mask_sine, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_1[:,i],RA_lines_2[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_1[:,i],dec_lines_2[:,i], color='white')
ax.set_xticks([       5,   35,   64,   93,  122,  150,  180])
ax.set_xticklabels(['60', '65', '70', '75', '80', '85', '95'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 22,  82, 142, 202])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sine_diff_false_B.pdf")


e1map_sterio, e2map_sterio, mask_sterio = DES.make_shear_rot_sterio_map(e1, e2_flip, c1, c2_flip, weight, mcorr, ra_flip, dec, 5.0, 2, Nx=250, Ny=250)

# plt.figure()
# imgplot = plt.imshow(e1map_sterio,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.title("e1_sterio")
# plt.figure()
# imgplot = plt.imshow(e2map_sterio,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.title("e2_sterio")

# plt.show()

gamma_plane = e1map_sterio + 1j*e2map_sterio

kappa_sterio = mm.reduced_shear_to_kappa_plane(gamma_plane, 1.0, 1.0, sigma=2*8./2.355*sigma_scale, Iterate=Iterate)

zoom_region = 2.0*np.arctan(125*5.0*np.pi/(60.*180.*2.0))
RA_lines_1, RA_lines_2, dec_lines_1, dec_lines_2 = sterio_grid_lines(resolution=250,zoom_region=zoom_region)

fig, ax = plt.subplots()
imgplot = ax.imshow(kappa_sterio.real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sterio, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sterio_E.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(-kappa_sterio.imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
# minus sign to correct for "mirror flip which effects B-mode"
ax.imshow(mask_sterio, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sterio_B.pdf")


# plot MW and setrio difference

fig, ax = plt.subplots()
imgplot = ax.imshow(kappa_sterio.real-k_mw_south_real,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
ax.imshow(mask_sterio, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.imshow(mask_south_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sterio_mw_diff_E.pdf")

fig, ax = plt.subplots()
imgplot = ax.imshow(-kappa_sterio.imag-k_mw_south_imag,interpolation='nearest',vmin=-0.015,vmax=0.015, cmap="cubehelix")
# minus sign to correct for "mirror flip which effects B-mode"
ax.imshow(mask_sterio, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.imshow(mask_south_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
ax.set_autoscale_on(False)
for i in range(RA_lines_2.shape[1]):
	ax.plot(RA_lines_2[:,i],RA_lines_1[:,i], color='white')
for i in range(dec_lines_2.shape[1]):
	ax.plot(dec_lines_2[:,i],dec_lines_1[:,i], color='white')
ax.set_xticks([      24,   53,   81,   108, 135,  163,  190,  219, 248])
ax.set_xticklabels(['55', '60', '65', '70', '75', '80', '85', '95', '100'])
ax.set_xlabel('RA (degrees)')
ax.set_yticks([ 41,  103, 165, 228])
ax.set_yticklabels(['-45', '-50', '-55', '-60'])
ax.set_ylabel('Dec (degrees)')
plt.colorbar(imgplot)
if save_figs:
	plt.savefig("fig/DES_sterio_mw_diff_B.pdf")

if show_figs:
	plt.show()
