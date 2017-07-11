#importing packages                                                                                                                                
import scipy.io as sio
import numpy as np
import scipy as sp
import scipy.ndimage as im
import random as rn
import healpy as hp
import matplotlib.mlab as mlab
import pyssht as ssht
import cy_mass_mapping as mm
import cy_healpy_mass_mapping as hp_mm
import cy_lsst_challenge_1_utils as lsst
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import cm
import os

do_hp = True
do_mw = False


N_side = 512
N_pix = 12*N_side*N_side
L_hp = N_side*2+1

short_cut_today_ready_hp = True

challenge_data_dir = "/Users/cwallis/Documents/Projects/LSST_mass_mapping/challenge_1_data/"

if do_hp:

	if not short_cut_today_ready_hp:
		print "Reading galaxy files..."
		filenames = [challenge_data_dir+"galaxy_files/"+f for f in os.listdir(challenge_data_dir+"galaxy_files/")]

		e1map, e2map, weight = lsst.make_shear_map_hp(filenames, N_side, N_mingal=10)

		np.save(challenge_data_dir+"maps/hp_e1_"+str(N_side)+".npy", e1map)
		np.save(challenge_data_dir+"maps/hp_e2_"+str(N_side)+".npy", e2map)
		print "done"
	else:
		print "Reading hp maps (N_side=" + str(N_side) + ")..."
		e1map = np.load(challenge_data_dir+"maps/hp_e1_"+str(N_side)+".npy")
		e2map = np.load(challenge_data_dir+"maps/hp_e2_"+str(N_side)+".npy")
		print "done"

	print "Smoothing shear maps"
	dummy, e1map_smooth, e2map_smooth = hp.sphtfunc.smoothing([np.zeros(N_pix),e1map,e2map], sigma=np.pi*0.5/180)

	print "Recovering kappa maps: HEALPix"
	kappa_E_hp, kappa_B_hp = hp_mm.gamma_to_kappa_hp(e1map_smooth, e2map_smooth, L_hp, N_side)

	indecies = np.where(e1map_smooth == hp.UNSEEN)

	kappa_E_hp[indecies] = hp.UNSEEN
	kappa_B_hp[indecies] = hp.UNSEEN

	# hp.mollview(e1map)
	# hp.mollview(e2map)
	hp.mollview(e1map_smooth, min=-0.01, max=0.01, cmap="cubehelix", title="e1 HEALPix")
	plt.savefig("fig/lsst_challenge_hp_e1.pdf")
	hp.mollview(e2map_smooth, min=-0.01, max=0.01, cmap="cubehelix", title="e2 HEALPix")
	plt.savefig("fig/lsst_challenge_hp_e2.pdf")
	hp.mollview(kappa_E_hp, min=-0.01, max=0.01, cmap="cubehelix", title="E HEALPix")
	plt.savefig("fig/lsst_challenge_hp_kappa_E.pdf")
	hp.mollview(kappa_B_hp, min=-0.01, max=0.01, cmap="cubehelix", title="B HEALPix")
	plt.savefig("fig/lsst_challenge_hp_kappa_B.pdf")
# plt.show()

if do_mw:

	L_mw = 1024

	short_cut_today_ready_mw = True

	challenge_data_dir = "/Users/cwallis/Documents/Projects/LSST_mass_mapping/challenge_1_data/"


	if not short_cut_today_ready_mw:
		print "Reading galaxy files..."
		filenames = [challenge_data_dir+"galaxy_files/"+f for f in os.listdir(challenge_data_dir+"galaxy_files/")]

		e1map, e2map, weight = lsst.make_shear_map_mw(filenames, L_mw, N_mingal=10)

		np.save(challenge_data_dir+"maps/mw_e1_"+str(L_mw)+".npy", e1map)
		np.save(challenge_data_dir+"maps/mw_e2_"+str(L_mw)+".npy", e2map)
		print "done"
	else:
		print "Reading mw maps (L_mw=" + str(L_mw) + ")..."
		e1map = np.load(challenge_data_dir+"maps/mw_e1_"+str(L_mw)+".npy")
		e2map = np.load(challenge_data_dir+"maps/mw_e2_"+str(L_mw)+".npy")
		print "done"


	emap = -e1map - 1j*e2map

	print "Recovering kappa maps: MW"
	kappa_mw = mm.gamma_to_kappa_mw(emap, L_mw, sigma=np.pi*0.5/180)

	indecies = np.where(emap != emap)
	emap[indecies] = 0.0
	elm = ssht.forward(emap, L_mw, Spin=2)
	elm_smooth = ssht.guassian_smoothing(elm, L_mw, sigma_in=np.pi*0.5/180)
	emap_smooth = -ssht.inverse(elm_smooth, L_mw, Spin=2)
	emap_smooth[indecies] = np.nan

	emap_smooth = np.fliplr(emap_smooth)
	emap_smooth = np.flipud(emap_smooth)
	f_plot_real, mask_array_real, f_plot_imag, mask_array_imag = ssht.mollweide_projection(emap_smooth, L_mw, resolution=1000)#, rot=[0.0, np.pi/4, np.pi/2])
	plt.figure()
	imgplot = plt.imshow(f_plot_real,interpolation='nearest', vmin=-0.01, vmax=0.01, cmap="cubehelix")
	plt.colorbar(imgplot,fraction=0.025, pad=0.04)
	plt.imshow(mask_array_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
	plt.gca().set_aspect("equal")
	plt.title("e1 MW")
	plt.axis('off')
	plt.savefig("fig/lsst_challenge_mw_e1.pdf")

	plt.figure()
	imgplot = plt.imshow(f_plot_imag,interpolation='nearest', vmin=-0.01, vmax=0.01, cmap="cubehelix")
	plt.colorbar(imgplot,fraction=0.025, pad=0.04)
	plt.imshow(mask_array_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
	plt.gca().set_aspect("equal")
	plt.title("e2 MW")
	plt.axis('off')
	plt.savefig("fig/lsst_challenge_mw_e2.pdf")

	kappa_mw = np.fliplr(kappa_mw)
	kappa_mw = np.flipud(kappa_mw)
	f_plot_real, mask_array_real, f_plot_imag, mask_array_imag = ssht.mollweide_projection(kappa_mw, L_mw, resolution=1000)#, rot=[0.0, np.pi/4, np.pi/2])
	plt.figure()
	imgplot = plt.imshow(f_plot_real,interpolation='nearest', vmin=-0.01, vmax=0.01, cmap="cubehelix")
	plt.colorbar(imgplot,fraction=0.025, pad=0.04)
	plt.imshow(mask_array_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
	plt.gca().set_aspect("equal")
	plt.title("E MW")
	plt.axis('off')
	plt.savefig("fig/lsst_challenge_mw_kappa_E.pdf")

	plt.figure()
	imgplot = plt.imshow(f_plot_imag,interpolation='nearest', vmin=-0.01, vmax=0.01, cmap="cubehelix")
	plt.colorbar(imgplot,fraction=0.025, pad=0.04)
	plt.imshow(mask_array_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
	plt.gca().set_aspect("equal")
	plt.title("B MW")
	plt.axis('off')
	plt.savefig("fig/lsst_challenge_mw_kappa_B.pdf")

plt.show()