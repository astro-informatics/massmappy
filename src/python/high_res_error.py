import numpy as np
import matplotlib.pyplot as plt
import pyssht as ssht
import massmappy.cy_mass_mapping as mm
from matplotlib import cm

def plot_projection(f, mask, title, save_figs, filename):
	plt.figure()
	imgplot = plt.imshow(f,interpolation='nearest')
	plt.colorbar(imgplot)
	plt.imshow(mask, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
	plt.title(title)
	plt.axis('off')
	if save_figs:
		plt.savefig(filename)

def error_function(x):
	error = np.sqrt(np.mean(x**2))
	return error

orth_resolution = 2000 # 2000 4000 for pixels similar size to nside 2048
L = orth_resolution*2
cyl_resolution = 2*L-1
Method = "MW"
save_figs = False #True
show_figs = False #True
do_cylindrical = True
do_full_sky_proj = True
do_half_sky_proj = True
Cyl_Projection_array = ["MERCATOR", "SINE"]
Projection_array  = ["OP",   "SP",         "GP"]
zoom_region_array = [np.pi/2-1E-6, np.pi/2, np.pi/4]
N_real = 10
#Projection_array = ["OP"]

print "calculate errors"
N_angle = 12*orth_resolution/500
angle, angle_step = np.linspace(0.0, np.pi/2, N_angle, retstep=True)
thetas, phis       = ssht.sample_positions(L, Method=Method)
x = np.linspace(-1., 1., orth_resolution)
y = np.linspace(-1., 1., orth_resolution)
x, y = np.meshgrid(x, y)
rho = np.sqrt(x*x + y*y)

error_cylindrical_E   = np.zeros((N_angle,N_real))
error_mercator_E      = np.zeros((N_angle,N_real))
error_sine_E          = np.zeros((N_angle,N_real))
error_orthographic_E  = np.zeros((N_angle,N_real))
error_steriographic_E = np.zeros((N_angle,N_real))
error_gnomic_E        = np.zeros((N_angle,N_real))

error_cylindrical_B   = np.zeros((N_angle,N_real))
error_mercator_B      = np.zeros((N_angle,N_real))
error_sine_B          = np.zeros((N_angle,N_real))
error_orthographic_B  = np.zeros((N_angle,N_real))
error_steriographic_B = np.zeros((N_angle,N_real))
error_gnomic_B        = np.zeros((N_angle,N_real))

Cls = np.loadtxt("data/cls_ap.txt")

for i_real in range(N_real):
	print i_real
	k_lm_mw = mm.generate_kappa_lm_mw(np.array(Cls[:,1]), L)

	ks_lm_mw = ssht.guassian_smoothing(k_lm_mw, L, sigma_in=np.pi/256)

	k_mw = ssht.inverse(ks_lm_mw, L, Reality=True, Method=Method)

	gamma_lm = mm.kappa_lm_to_gamma_lm_mw(ks_lm_mw, L)

	gamma = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)
	if do_cylindrical:
		print "Doing Projection Cylindrical"

		gamma_plane = -gamma
		kappa_orig_plane = k_mw
		kappa_plane  = mm.gamma_to_kappa_plane(gamma_plane, np.pi/L, 2*np.pi/(2*L-1))

		for i in range(N_angle):
			x = np.abs(thetas-np.pi/2)
			indexes = np.nonzero((x < angle[i]+angle_step/2) & (x > angle[i]-angle_step/2))
			error_cylindrical_E[i,i_real] = error_function((k_mw[indexes]-kappa_plane[indexes]).real)/error_function(k_mw[indexes])
			error_cylindrical_B[i,i_real] = error_function((kappa_plane[indexes]).imag)/error_function(k_mw[indexes])

	if do_full_sky_proj:

		for Projection in Cyl_Projection_array:

			print "Doing Projection " + Projection

			proj_real, mask_real, proj_imag, mask_imag,\
				 = ssht.equatorial_projection(gamma, L, resolution=cyl_resolution, Method=Method, Projection=Projection, \
				 	Spin=2)

			gamma_plane = - proj_real + 1j*proj_imag

			# run planar KS
			kappa_plane  = mm.gamma_to_kappa_plane(gamma_plane, 1.0,1.0)

			kappa_orig, mask \
				= ssht.equatorial_projection(k_mw, L, resolution=cyl_resolution, Method=Method, Projection=Projection)

			theta_angle, phi_angle = ssht.equatorial_projection_angle_array(cyl_resolution, Projection=Projection)
			phi_angle[phi_angle>np.pi] = phi_angle[phi_angle>np.pi]-2*np.pi

			for i in range(N_angle):
				if Projection == "MERCATOR":
					x = np.abs(theta_angle-np.pi/2)
					y = np.abs(phi_angle)
					indexes = np.nonzero((x < angle[i]+angle_step/2) & (x > angle[i]-angle_step/2))# & (y < angle[i]+angle_step) & (y > angle[i]-angle_step))
					error_mercator_E[i,i_real] = error_function((kappa_orig[indexes]-kappa_plane[indexes]).real)/error_function(kappa_orig[indexes])
					error_mercator_B[i,i_real] = error_function((kappa_plane[indexes]).imag)/error_function(kappa_orig[indexes])
				if Projection == "SINE":
					x = np.sqrt((theta_angle-np.pi/2)*(theta_angle-np.pi/2)+(phi_angle)*(phi_angle))
					indexes = np.nonzero((x < angle[i]+angle_step/2) & (x > angle[i]-angle_step/2))
					error_sine_E[i,i_real] = error_function((kappa_orig[indexes]-kappa_plane[indexes]).real)/error_function(kappa_orig[indexes])
					error_sine_B[i,i_real] = error_function((kappa_plane[indexes]).imag)/error_function(kappa_orig[indexes])
			# 	index_region = np.zeros(x.shape)
			# 	index_region[indexes] = 1.0
			# 	plt.figure()
			# 	plt.imshow(index_region)
			# 	plt.show()
			# plt.figure()
			# im = plt.imshow(kappa_orig- kappa_plane.real)
			# plt.colorbar(im)
			# plt.figure()
			# im = plt.imshow(kappa_plane.imag)
			# plt.colorbar(im)
			# plt.show()

	if do_half_sky_proj:

		for Projection, zoom_region in zip(Projection_array, zoom_region_array):
			print "Doing Projection  ", Projection

			# project gamma
			proj_north_real, mask_north_real, proj_south_real, mask_south_real,\
			proj_north_imag, mask_north_imag, proj_south_imag, mask_south_imag\
				 = ssht.polar_projection(gamma, L, resolution=orth_resolution, Method=Method, rot=[0.0,np.pi/2,0.0], \
				 	Projection=Projection, Spin=2, zoom_region=zoom_region)

			gamma_plane_north = proj_north_real + 1j*proj_north_imag


			# Project original
			kappa_orig_north, mask_north, kappa_orig_south, mask_south \
				= ssht.polar_projection(k_mw, L, resolution=orth_resolution, Method=Method, rot=[0.0,np.pi/2,0.0], \
					Projection=Projection, zoom_region=zoom_region)

			# run planar KS
			kappa_plane_north  = mm.gamma_to_kappa_plane(gamma_plane_north, 1.0, 1.0)



			for i in range(N_angle):
				if Projection == "OP":
					indexes = np.nonzero((rho < np.sin(min(angle[i]+angle_step/2,np.pi/2))) &\
								(rho > np.sin(angle[i]-angle_step/2)))
					error_orthographic_E[i,i_real] = error_function(kappa_orig_north[indexes]-kappa_plane_north[indexes].real)/error_function(kappa_orig_north[indexes])
					error_orthographic_B[i,i_real] = error_function(kappa_plane_north[indexes].imag)/error_function(kappa_orig_north[indexes])
				if Projection == "SP":
					indexes = np.nonzero((rho < np.tan(min(angle[i]+angle_step/2,np.pi/2)/2.0)) &\
								(rho > np.tan((angle[i]-angle_step/2)/2.0)))
					error_steriographic_E[i,i_real] = error_function(kappa_orig_north[indexes]-kappa_plane_north[indexes].real)/error_function(kappa_orig_north[indexes])
					error_steriographic_B[i,i_real] = error_function(kappa_plane_north[indexes].imag)/error_function(kappa_orig_north[indexes])
				if Projection == "GP":
					if angle[i] < np.pi/4:
						indexes = np.nonzero((rho < np.tan(min(angle[i]+angle_step/2,np.pi/4))) &\
									(rho > np.sin(angle[i]-angle_step/2)))
						error_gnomic_E[i,i_real] = error_function(kappa_orig_north[indexes]-kappa_plane_north[indexes].real)/error_function(kappa_orig_north[indexes])
						error_gnomic_B[i,i_real] = error_function(kappa_plane_north[indexes].imag)/error_function(kappa_orig_north[indexes])
					else:
						error_gnomic_E[i,i_real] = np.nan
						error_gnomic_B[i,i_real] = np.nan



error_cylindrical_E_av   = error_cylindrical_E.mean(axis=1)
error_mercator_E_av      = error_mercator_E.mean(axis=1)
error_sine_E_av          = error_sine_E.mean(axis=1)
error_orthographic_E_av  = error_orthographic_E.mean(axis=1)
error_steriographic_E_av = error_steriographic_E.mean(axis=1)
error_gnomic_E_av        = error_gnomic_E.mean(axis=1)

error_cylindrical_E_std   = error_cylindrical_E.std(axis=1)
error_mercator_E_std      = error_mercator_E.std(axis=1)
error_sine_E_std          = error_sine_E.std(axis=1)
error_orthographic_E_std  = error_orthographic_E.std(axis=1)
error_steriographic_E_std = error_steriographic_E.std(axis=1)
error_gnomic_E_std        = error_gnomic_E.std(axis=1)

error_cylindrical_B_av   = error_cylindrical_B.mean(axis=1)
error_mercator_B_av      = error_mercator_B.mean(axis=1)
error_sine_B_av          = error_sine_B.mean(axis=1)
error_orthographic_B_av  = error_orthographic_B.mean(axis=1)
error_steriographic_B_av = error_steriographic_B.mean(axis=1)
error_gnomic_B_av        = error_gnomic_B.mean(axis=1)

error_cylindrical_B_std   = error_cylindrical_B.std(axis=1)
error_mercator_B_std      = error_mercator_B.std(axis=1)
error_sine_B_std          = error_sine_B.std(axis=1)
error_orthographic_B_std  = error_orthographic_B.std(axis=1)
error_steriographic_B_std = error_steriographic_B.std(axis=1)
error_gnomic_B_std        = error_gnomic_B.std(axis=1)

if do_cylindrical:
	np.savetxt("data/high_res_results/error_cylindrical_E.txt", error_cylindrical_E)
if do_full_sky_proj:
	np.savetxt("data/high_res_results/error_mercater_E.txt", error_mercator_E)
	np.savetxt("data/high_res_results/error_sine_E.txt", error_sine_E)
if do_half_sky_proj:
	np.savetxt("data/high_res_results/error_orthographic_E.txt", error_orthographic_E)
	np.savetxt("data/high_res_results/error_steriographic_E.txt", error_steriographic_E)
	np.savetxt("data/high_res_results/error_gnomic_E.txt", error_gnomic_E)

if do_cylindrical:
	np.savetxt("data/high_res_results/error_cylindrical_B.txt", error_cylindrical_B)
if do_full_sky_proj:
	np.savetxt("data/high_res_results/error_mercater_B.txt", error_mercator_B)
	np.savetxt("data/high_res_results/error_sine_B.txt", error_sine_B)
if do_half_sky_proj:
	np.savetxt("data/high_res_results/error_orthographic_B.txt", error_orthographic_B)
	np.savetxt("data/high_res_results/error_steriographic_B.txt", error_steriographic_B)
	np.savetxt("data/high_res_results/error_gnomic_B.txt", error_gnomic_B)


CFHTLens = np.sqrt(154.)
DES_SV = np.sqrt(139.)
KiDS = np.sqrt(1500.)
DES_full = np.sqrt(5000.)
Euclid = np.sqrt(15000.)

label_buffer = 2.0
ell = 200
#normalisation = np.sqrt(Cls[ell,1]*ell*(ell+1)/(2*np.pi))
normalisation = 1.0

plt.figure()
plt.plot([CFHTLens,CFHTLens], [0,1], 'k--')
plt.text(CFHTLens+label_buffer, 0.55, "CFHTLenS", fontsize=18)
plt.plot([DES_full,DES_full], [0,1], 'k--')
plt.text(DES_full+label_buffer, 0.9, "DES full", fontsize=18)
plt.plot([Euclid, Euclid], [0,1], 'k--')
plt.text(Euclid+label_buffer, 0.9, "Euclid", fontsize=18)
plt.errorbar(angle*2*180/np.pi, error_cylindrical_E_av/normalisation,color='red',\
	yerr=error_cylindrical_E_std/normalisation, label="Cylindrical")
plt.errorbar(angle*2*180/np.pi, error_mercator_E_av/normalisation,color='magenta',\
	yerr=error_mercator_E_std/normalisation, label="Mercator")
plt.errorbar(angle*2*180/np.pi, error_sine_E_av/normalisation,color='cyan',\
	yerr=error_sine_E_std/normalisation, label="Sinusiodal")
plt.errorbar(angle*2*180/np.pi, error_orthographic_E_av/normalisation,color='green',\
	yerr=error_orthographic_E_std/normalisation, label="Orthographic")
plt.errorbar(angle*2*180/np.pi, error_steriographic_E_av/normalisation,color='blue',\
	yerr=error_steriographic_E_std/normalisation, label="Steriographic")
plt.errorbar(angle*2*180/np.pi, error_gnomic_E_av/normalisation,color='brown',\
	yerr=error_gnomic_E_std/normalisation, label="Gnomic")
plt.xlabel(r"$\theta  ({\rm degrees})$", fontsize=20)
plt.ylabel(r"$\sigma_{\rm err}/\sigma_{\kappa}$", fontsize=20)
plt.legend(loc=2)
plt.xlim([0,180])
plt.ylim([0,1.0])
if save_figs:
	plt.savefig("fig/ring_error_E.png")

plt.figure()
plt.plot([CFHTLens,CFHTLens], [0,1], 'k--')
plt.text(CFHTLens+label_buffer, 0.55, "CFHTLenS", fontsize=18)
plt.plot([DES_full,DES_full], [0,1], 'k--')
plt.text(DES_full+label_buffer, 0.9, "DES full", fontsize=18)
plt.plot([Euclid, Euclid], [0,1], 'k--')
plt.text(Euclid+label_buffer, 0.9, "Euclid", fontsize=18)
plt.errorbar(angle*2*180/np.pi, error_cylindrical_B_av/normalisation,color='r',\
	yerr=error_cylindrical_B_std/normalisation, label="Cylindrical")
plt.errorbar(angle*2*180/np.pi, error_mercator_B_av/normalisation,color='magenta',\
	yerr=error_mercator_B_std/normalisation, label="Mercator")
plt.errorbar(angle*2*180/np.pi, error_sine_B_av/normalisation,color='cyan',\
	yerr=error_sine_B_std/normalisation, label="Sinusiodal")
plt.errorbar(angle*2*180/np.pi, error_orthographic_B_av/normalisation,color='green',\
	yerr=error_orthographic_B_std/normalisation, label="Orthographic")
plt.errorbar(angle*2*180/np.pi, error_steriographic_B_av/normalisation,color='blue',\
	yerr=error_steriographic_B_std/normalisation, label="Steriographic")
plt.errorbar(angle*2*180/np.pi, error_gnomic_B_av/normalisation,color='brown',\
	yerr=error_gnomic_B_std/normalisation, label="Gnomic")
plt.xlabel(r"$\theta  ({\rm degrees})$", fontsize=20)
plt.ylabel(r"$\sigma_{\rm err}/\sigma_{\kappa}$", fontsize=20)
plt.legend(loc=2)
plt.xlim([0,180])
plt.ylim([0,1.0])
if save_figs:
	plt.savefig("fig/ring_error_B.png")

if show_figs:
	plt.show()
