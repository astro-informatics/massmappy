import numpy as np
import scipy.ndimage as spim
import matplotlib.pyplot as plt
import pyssht as ssht
import massmappy.cy_mass_mapping as mm
import massmappy.cy_healpy_mass_mapping as hp_mm
import healpy as hp
from matplotlib import cm
from matplotlib import ticker
import time

def plot_projection(f, mask, title, plot_titles, save_figs, filename, smooth=-1, color_range=0.016):
	plt.figure()
	if smooth>-0.5:
		f = spim.gaussian_filter(f, smooth)
	imgplot = plt.imshow(f,interpolation='nearest',vmin=-color_range,vmax=color_range, cmap="cubehelix")
	cb = plt.colorbar(imgplot, orientation='horizontal')
	tick_locator = ticker.MaxNLocator(nbins=8)
	cb.locator = tick_locator
	cb.update_ticks()
	plt.imshow(mask, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
	if plot_titles:
		plt.title(title)
	plt.axis('off')
	if save_figs:
		plt.savefig(filename)


L = 512
orth_resolution = 500
cyl_resolotion = 2*L-1
mol_resolution = 500
Method = "MW"
save_figs = True
show_figs = False
plot_titles = True
Equatorial_Projection_array = ["MERCATOR", "SINE"]
#Cyl_Projection_array = ["SINE"]
Polar_Projection_array = ["OP", "SP", "GP"]
#Projection_array = ["SP"]
smooth = -1
kappa_plot_range = 0.016
gamma_plot_range = 0.01

Iterate  = True
tol_error = 1E-10
print_iteration_results = True

Cls = np.loadtxt("data/cls_ap.txt")

print "mm.generate_kappa_lm_mw"
k_lm_mw = mm.generate_kappa_lm_mw(np.array(Cls[:,1]), L, 2)

print "ssht.guassian_smoothing"
ks_lm_mw = ssht.guassian_smoothing(k_lm_mw, L, sigma_in=np.pi/256)

print "ssht.inverse"
k_mw = ssht.inverse(ks_lm_mw, L, Reality=True, Method=Method)

print "mm.kappa_lm_to_gamma_lm_mw"
gamma_lm = mm.kappa_lm_to_gamma_lm_mw(ks_lm_mw, L)

print "ssht.inverse"
gamma = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)

print "making shear"
if Iterate:
	shear = gamma/(1-k_mw)
else:
	shear = gamma.copy()

if print_iteration_results:
	ell = 200
	normalisation = np.sqrt(Cls[ell,1]*ell*(ell+1)/(2*np.pi))
	f = open('data/reduced_shear_iteration.txt','w')
	print "getting info on spherical MW inversion"
	start_mw = time.clock()
	k_mw_rec, count_mw = mm.reduced_shear_to_kappa_mw(shear, L, Method=Method, tol_error=tol_error, Iterate=Iterate, return_count=True)
	elapsed_mw = (time.clock() - start_mw)

	error_mw = np.sqrt(np.mean((k_mw-k_mw_rec.real)**2))
	f.write('Spherical MW ' + str(count_mw) + '  ' + str(error_mw/normalisation) + '  ' + str(elapsed_mw) + '\n')

	print "getting info on spherical MW inversion"
	print "mm.lm_hp2lm"
	ks_lm_hp = mm.lm2lm_hp(ks_lm_mw, L)
	Nside = L/2
	print "hp.alm2map"
	k_hp = hp.alm2map(ks_lm_hp,Nside, lmax=L-1)
	print "mm.kappa_lm_to_gamma_lm_hp"
	gamma_E_lm_hp, gamma_B_lm_hp = mm.kappa_lm_to_gamma_lm_hp(ks_lm_hp,L)
	print "hp.alm2map"
	gamma_lm_list = [ks_lm_hp, gamma_E_lm_hp, gamma_B_lm_hp]
	maps_hp = hp.alm2map(gamma_lm_list, Nside)
	start_hp = time.clock()
	kappa_E_map_hp_rec, kappa_B_map_hp_rec, count_hp = hp_mm.reduced_shear_to_kappa_hp(maps_hp[1], maps_hp[2], L, Nside, tol_error=tol_error, Iterate=Iterate, return_count=True)
	elapsed_hp = (time.clock() - start_hp)

	error_hp = np.sqrt(np.mean((k_hp-kappa_E_map_hp_rec)**2))
	f.write('Spherical healpy ' + str(count_hp) + '  ' + str(error_hp/normalisation) + '  ' + str(elapsed_hp) + '\n')

for Projection in Equatorial_Projection_array:

	print "Doing Projection " + Projection

	proj_real, mask_real, proj_imag, mask_imag,\
		 = ssht.equatorial_projection(shear, L, resolution=cyl_resolotion, Method=Method, Projection=Projection, \
		 	Spin=2)

	shear_plane = - proj_real + 1j*proj_imag


	# Project original
	kappa_orig, mask \
		= ssht.equatorial_projection(k_mw, L, resolution=cyl_resolotion, Method=Method, Projection=Projection)

	if print_iteration_results:
		print "getting info on planar inversion"
		start_plane = time.clock()
		kappa_plane, count_plane  = mm.reduced_shear_to_kappa_plane(shear_plane, 1.0,1.0, tol_error=tol_error, Iterate=Iterate, return_count=True)
		elapsed_plane = (time.clock() - start_plane)

		error_plane = np.sqrt(np.nanmean((kappa_orig-kappa_plane.real)**2))
		print count_plane, error_plane, elapsed_plane
		f.write(Projection + ' ' + str(count_plane) + '  ' + str(error_plane/normalisation) + '  ' + str(elapsed_plane) + '\n')
	else:
		# run planar KS
		kappa_plane  = mm.reduced_shear_to_kappa_plane(shear_plane, 1.0,1.0, tol_error=tol_error, Iterate=Iterate)

		plot_projection(kappa_orig, mask, "Input kappa plane north "+Projection,\
			plot_titles, save_figs, "fig/low_res_input_kappa_"+Projection+".pdf", color_range=kappa_plot_range)

		plot_projection(kappa_plane.real, mask, "kappa plane real "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_real.pdf", color_range=kappa_plot_range)
		plot_projection(kappa_plane.imag, mask, "kappa plane imag "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_imag.pdf", color_range=kappa_plot_range)

		plot_projection(kappa_plane.real-kappa_orig, mask, "kappa plane error real "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_error_real.pdf", smooth=smooth, color_range=kappa_plot_range)
		plot_projection(kappa_plane.imag, mask, "kappa plane error imag "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_error_imag.pdf", smooth=smooth, color_range=kappa_plot_range)

		plot_projection(shear_plane.real, mask_real, "gamma plane real "+Projection,\
			plot_titles, save_figs, "fig/low_res_gamma_"+Projection+"_real.pdf", color_range=gamma_plot_range)
		plot_projection(shear_plane.imag, mask_imag, "gamma plane imag "+Projection,\
			plot_titles, save_figs, "fig/low_res_gamma_"+Projection+"_imag.pdf", color_range=gamma_plot_range)

if show_figs:
	plt.show()

print "Doing Projection Simple Cylindrical"

shear_plane = -shear
kappa_orig_plane = k_mw

if print_iteration_results:
	print "getting info on planar inversion"
	start_plane = time.clock()
	kappa_plane, count_plane  = mm.reduced_shear_to_kappa_plane(shear_plane, 1.0,1.0, tol_error=tol_error, Iterate=Iterate, return_count=True)
	elapsed_plane = (time.clock() - start_plane)

	error_plane = np.sqrt(np.nanmean((kappa_orig_plane-kappa_plane.real)**2))
	print count_plane, error_plane, elapsed_plane
	f.write('Cylindrical ' + str(count_plane) + '  ' + str(error_plane/normalisation) + '  ' + str(elapsed_plane) + '\n')
else:
	kappa_plane  = mm.reduced_shear_to_kappa_plane(shear_plane, np.pi/L, 2*np.pi/(2*L-1), tol_error=tol_error, Iterate=Iterate)


	plt.figure()
	imgplot = plt.imshow((kappa_orig_plane),interpolation='nearest', extent=[0,2*np.pi,0,np.pi],\
		vmin=-kappa_plot_range,vmax=kappa_plot_range, cmap="cubehelix")
	plt.colorbar(imgplot, orientation='horizontal')
	if plot_titles:
		plt.title("Input kappa plane Cylindrical")
	plt.xlabel(r"$\phi$")
	plt.ylabel(r"$\theta$")
	if save_figs:
		plt.savefig("fig/low_res_input_kappa_cylindrical.pdf")

	plt.figure()
	imgplot = plt.imshow(kappa_plane.real,interpolation='nearest', extent=[0,2*np.pi,0,np.pi],\
		vmin=-kappa_plot_range,vmax=kappa_plot_range, cmap="cubehelix")
	plt.colorbar(imgplot, orientation='horizontal')
	if plot_titles:
		plt.title("kappa plane north real Cylindrical")
	plt.xlabel(r"$\phi$")
	plt.ylabel(r"$\theta$")
	if save_figs:
		plt.savefig("fig/low_res_kappa_cylindrical_real.pdf")


	plt.figure()
	imgplot = plt.imshow((kappa_orig_plane-kappa_plane).real,interpolation='nearest', extent=[0,2*np.pi,0,np.pi],\
		vmin=-kappa_plot_range,vmax=kappa_plot_range, cmap="cubehelix")
	plt.colorbar(imgplot, orientation='horizontal')
	if plot_titles:
		plt.title("kappa plane north error real Cylindrical")
	plt.xlabel(r"$\phi$")
	plt.ylabel(r"$\theta$")
	if save_figs:
		plt.savefig("fig/low_res_kappa_Cylindrical_north_error_real.pdf")


	plt.figure()
	imgplot = plt.imshow((kappa_orig_plane-kappa_plane).imag,interpolation='nearest', extent=[0,2*np.pi,0,np.pi],\
		vmin=-kappa_plot_range,vmax=kappa_plot_range, cmap="cubehelix")
	plt.colorbar(imgplot, orientation='horizontal')
	if plot_titles:
		plt.title("kappa plane north error imag Cylindrical")
	plt.xlabel(r"$\phi$")
	plt.ylabel(r"$\theta$")
	if save_figs:
		plt.savefig("fig/low_res_kappa_Cylindrical_north_error_imag.pdf")


	plt.figure()
	imgplot = plt.imshow(shear_plane.imag,interpolation='nearest', extent=[0,2*np.pi,0,np.pi],\
		vmin=-gamma_plot_range,vmax=gamma_plot_range, cmap="cubehelix")
	plt.colorbar(imgplot, orientation='horizontal')
	if plot_titles:
		plt.title("gamma plane imag Cylindrical")
	plt.xlabel(r"$\phi$")
	plt.ylabel(r"$\theta$")
	if save_figs:
		plt.savefig("fig/low_res_gamma_Cylindrical_imag.pdf")


	plt.figure()
	imgplot = plt.imshow(shear_plane.real,interpolation='nearest', extent=[0,2*np.pi,0,np.pi],\
		vmin=-gamma_plot_range,vmax=gamma_plot_range, cmap="cubehelix")
	plt.colorbar(imgplot, orientation='horizontal')
	if plot_titles:
		plt.title("gamma plane real Cylindrical")
	plt.xlabel(r"$\phi$")
	plt.ylabel(r"$\theta$")
	if save_figs:
		plt.savefig("fig/low_res_gamma_Cylindrical_real.pdf")

	if show_figs:
		plt.show()

plt.close("all")

for Projection in Polar_Projection_array:
	print "Doing Projection  ", Projection

	# project gamma
	proj_north_real, mask_north_real, proj_south_real, mask_south_real,\
	proj_north_imag, mask_north_imag, proj_south_imag, mask_south_imag\
		 = ssht.polar_projection(shear, L, resolution=orth_resolution, Method=Method, rot=[0.0,np.pi/2,0.0], Projection=Projection, Spin=2)

	shear_plane_north = proj_north_real + 1j*proj_north_imag
	shear_plane_south = proj_south_real - 1j*proj_south_imag

	# Project original
	kappa_orig_north, mask_north, kappa_orig_south, mask_south \
		= ssht.polar_projection(k_mw, L, resolution=orth_resolution, Method=Method, rot=[0.0,np.pi/2,0.0], Projection=Projection)

	if print_iteration_results:
		print "getting info on planar inversion"
		start_plane = time.clock()
		kappa_plane_north, count_plane  = mm.reduced_shear_to_kappa_plane(shear_plane_north, 1.0,1.0, tol_error=tol_error, Iterate=Iterate, return_count=True)
		elapsed_plane = (time.clock() - start_plane)

		error_plane = np.sqrt(np.nanmean((kappa_orig_north-kappa_plane_north.real)**2))
		print count_plane, error_plane, elapsed_plane
		f.write(Projection + ' ' + str(count_plane) + '  ' + str(error_plane/normalisation) + '  ' + str(elapsed_plane) + '\n')
	else:
		# run planar KS
		kappa_plane_north  = mm.reduced_shear_to_kappa_plane(shear_plane_north, 1.0, 1.0, tol_error=tol_error, Iterate=Iterate)
		kappa_plane_south  = mm.reduced_shear_to_kappa_plane(shear_plane_south, 1.0, 1.0, tol_error=tol_error, Iterate=Iterate)

		plot_projection(kappa_orig_north, mask_north_real, "Input kappa plane north "+Projection,\
			plot_titles, save_figs, "fig/low_res_input_kappa_"+Projection+"_north.pdf", color_range=kappa_plot_range)
		plot_projection(kappa_orig_south.real, mask_north_real, "Input kappa plane south "+Projection, \
			plot_titles, save_figs, "fig/low_res_input_kappa_"+Projection+"_south.pdf", color_range=kappa_plot_range)

		plot_projection(kappa_plane_north.real, mask_north_real, "kappa plane north real "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_north_real.pdf", color_range=kappa_plot_range)

		plot_projection(kappa_plane_north.real-kappa_orig_north, mask_north_real, "kappa plane north error real "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_north_error_real.pdf", smooth=smooth, color_range=kappa_plot_range)
		plot_projection(kappa_plane_north.imag, mask_north_imag, "kappa plane north error imag "+Projection,\
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_north_error_imag.pdf", smooth=smooth, color_range=kappa_plot_range)

		plot_projection(kappa_plane_south.real, mask_north_real, "kappa plane south real "+Projection, \
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_south_real.pdf", color_range=kappa_plot_range)

		plot_projection(kappa_plane_south.real-kappa_orig_south, mask_north_real, "kappa plane south error real "+Projection, \
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_south_error_real.pdf", smooth=smooth, color_range=kappa_plot_range)
		plot_projection(kappa_plane_south.imag, mask_north_imag, "kappa plane south error imag "+Projection, \
			plot_titles, save_figs, "fig/low_res_kappa_"+Projection+"_south_error_imag.pdf", smooth=smooth, color_range=kappa_plot_range)

		plot_projection(shear_plane_north.real, mask_north_real, "gamma plane north real "+Projection,\
			plot_titles, save_figs, "fig/low_res_gamma_"+Projection+"_north_real.pdf", color_range=gamma_plot_range)
		plot_projection(shear_plane_north.imag, mask_north_imag, "gamma plane north imag "+Projection,\
			plot_titles, save_figs, "fig/low_res_gamma_"+Projection+"_north_imag.pdf", color_range=gamma_plot_range)

		plot_projection(shear_plane_south.real, mask_north_real, "gamma plane south real "+Projection, \
			plot_titles, save_figs, "fig/low_res_gamma_"+Projection+"_south_real.pdf", color_range=gamma_plot_range)
		plot_projection(shear_plane_south.imag, mask_north_imag, "gamma plane south imag "+Projection, \
			plot_titles, save_figs, "fig/low_res_gamma_"+Projection+"_south_imag.pdf", color_range=gamma_plot_range)

		rotation_angle_north, rotation_angle_south = ssht.polar_projection_rotation_array(orth_resolution, Projection=Projection, rot=[0.0,np.pi/2,0.0])
		plot_projection(rotation_angle_north, mask_north_real, "rotation angles north "+Projection, \
			False, False, "dummy")
		plot_projection(rotation_angle_south, mask_south_real, "rotation angles south "+Projection, \
			False, False, "dummy")


	if show_figs:
		plt.show()
	plt.close("all")
