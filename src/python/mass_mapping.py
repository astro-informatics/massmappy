import numpy as np
import matplotlib.pyplot as plt
import pyssht as ssht
import cy_mass_mapping as mm
import cy_healpy_mass_mapping as hp_mm
from matplotlib import cm, colors, colorbar, gridspec
import healpy as hp


L = 257
lmax = L-1
orth_resolution = 500
Method = "MW"
Nside = 128
save_figs = True
show_figs = False

Cls = np.loadtxt("data/cls_ap.txt")

print "mm.generate_kappa_lm_mw"
k_lm_mw = mm.generate_kappa_lm_mw(np.array(Cls[:,1]), L, 2)

print "ssht.guassian_smoothing"
ks_lm_mw = ssht.guassian_smoothing(k_lm_mw, L, sigma_in=np.pi/L)

print "mm.lm_hp2lm"
ks_lm_hp = mm.lm2lm_hp(ks_lm_mw, L)

print "ssht.inverse"
k_mw = ssht.inverse(ks_lm_mw, L, Reality=True, Method=Method)

print "hp.alm2map"
k_hp = hp.alm2map(ks_lm_hp,Nside, lmax=lmax)

print "mm.kappa_lm_to_gamma_lm_mw"
gamma_lm = mm.kappa_lm_to_gamma_lm_mw(ks_lm_mw, L)

print "mm.kappa_lm_to_gamma_lm_hp"
gamma_E_lm_hp, gamma_B_lm_hp = mm.kappa_lm_to_gamma_lm_hp(ks_lm_hp,L)

print "ssht.inverse"
gamma = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)

print "hp.alm2map"
gamma_lm_list = [ks_lm_hp, gamma_E_lm_hp, gamma_B_lm_hp]
maps_hp = hp.alm2map(gamma_lm_list, Nside)

mask = np.zeros((L,2*L-1), dtype=float)#+1.0
mask[3*L/8:5*L/8,3*L/4:5*L/4] = 1.0
#gamma = gamma*mask

B = 15
mask_small = np.zeros((L,2*L-1), dtype=float)#+1.0
mask_small[3*L/8+B:5*L/8-B,3*L/4+B:5*L/4-B] = 1.0

print "Reconstruct on the sphere MW"
k_rec_mw = mm.gamma_to_kappa_mw(gamma, L, Method=Method)

print "Reconstruct on the sphere hp"

kappa_map_hp_rec_boris = hp_mm.gamma_to_kappa_hp_boris(maps_hp[1], maps_hp[2], L, Nside)
kappa_E_map_hp_rec, kappa_B_map_hp_rec = hp_mm.gamma_to_kappa_hp(maps_hp[1], maps_hp[2], L, Nside)

hp.mollview(k_hp, title="Kappa Healpis Input")
if save_figs:
	plt.savefig("fig/kappa_input_hp.png")
hp.mollview(maps_hp[1], title="Gamma Healpix real")
if save_figs:
	plt.savefig("fig/gamma_real_input_hp.png")

hp.mollview(maps_hp[2], title="Gamma Healpix imag")
if save_figs:
	plt.savefig("fig/gamma_imag_input_hp.png")

hp.mollview(kappa_E_map_hp_rec, title="Kappa Healpis Recovered")
if save_figs:
	plt.savefig("fig/kappa_rec_hp.png")

hp.mollview(kappa_E_map_hp_rec-k_hp, title="Kappa Healpis Recovered error real")
if save_figs:
	plt.savefig("fig/kappa_rec_error_real_hp.png")

hp.mollview(kappa_B_map_hp_rec, title="Kappa Healpis Recovered error imag")
if save_figs:
	plt.savefig("fig/kappa_rec_error_imag_hp.png")

hp.mollview(kappa_map_hp_rec_boris-k_hp, title="Kappa Healpis Recovered error boris chris")
if save_figs:
	plt.savefig("fig/kappa_rec_error_real_boris_hp.png")


print "on the plane cylindrical projection"

gamma_plane = -gamma#*mask
kappa_orig_plane = k_mw#*mask
kappa_plane  = mm.gamma_to_kappa_plane(gamma_plane, np.pi/L, 2*np.pi/(2*L-1))

print "on the plane orthographic projection"

orth_proj_north_real, mask_north_real, orth_proj_south_real, mask_south_real,\
orth_proj_north_imag, mask_north_imag, orth_proj_south_imag, mask_south_imag\
		 = ssht.orthographic_projection(gamma, L, resolution=orth_resolution, Method=Method, rot=[0.0,np.pi/2,0.0])
gamma_plane_orth_north = orth_proj_north_real + 1j*orth_proj_north_imag
gamma_plane_orth_south = orth_proj_south_real - 1j*orth_proj_south_imag
kappa_orig_orth_north, mask_north_real, kappa_orig_orth_south, mask_south_real \
	= ssht.orthographic_projection(k_mw, L, resolution=orth_resolution, Method=Method, rot=[0.0,np.pi/2,0.0])# make projection
kappa_plane_orth_north  = mm.gamma_to_kappa_plane(gamma_plane_orth_north, 1.0, 1.0)
kappa_plane_orth_south  = mm.gamma_to_kappa_plane(gamma_plane_orth_south, 1.0, 1.0)

# calculate errors
print "calculate errors"
N_angle = 12
angle, angle_step = np.linspace(0.0, np.pi/2, N_angle, retstep=True)
thetas, phis       = ssht.sample_positions(L, Method=Method)
x = np.linspace(-1., 1., orth_resolution)
y = np.linspace(-1., 1., orth_resolution)
x, y = np.meshgrid(x, y)
rho = np.sqrt(x*x + y*y)
error_sphere       = np.empty(N_angle)
error_cylindrical  = np.empty(N_angle)
error_orthographic = np.empty(N_angle)

for i in range(N_angle):
	indexes = np.abs(thetas-np.pi/2-angle[i])<angle_step/2
	error_sphere[i] = np.std(k_mw[indexes]-k_rec_mw[indexes])
	error_cylindrical[i] = np.std(k_mw[indexes]-kappa_plane[indexes])

#	indexes = rho-np.sin(angle[i]) < np.sin(angle[i]+angle_step) - np.sin(angle[i])
#	indexes = indexes and rho-np.sin(angle[i]) > np.sin(angle[i]-angle_step) - np.sin(angle[i])
	indexes = np.nonzero((rho < np.sin(min(angle[i]+angle_step/2,np.pi/2))) &\
				(rho > np.sin(angle[i]-angle_step/2)))
	error_orthographic[i] = np.std(kappa_orig_orth_north[indexes]-kappa_plane_orth_north[indexes])
	# print np.sin(angle[i]-angle_step/2), np.sin(min(angle[i]+angle_step/2,np.pi/2))
	# print rho[indexes].min(), rho[indexes].max()
	# print ""
	# kappa_plane_orth_south[:]=0.0
	# kappa_plane_orth_south[indexes]=1.0
	# plt.figure()
	# imgplot = plt.imshow(kappa_plane_orth_south.real,interpolation='nearest')#kappa_plane_orth_south.real
	# plt.colorbar(imgplot)
	# plt.axis('off')
	# plt.show()


plt.figure()
plt.plot(angle, error_sphere/np.std(k_mw))
plt.plot(angle, error_cylindrical/np.std(k_mw))
plt.plot(angle, error_orthographic/np.std(k_mw))
plt.xlabel(r"$\theta  ({\rm radians})$", fontsize=20)
plt.ylabel(r"$\sigma_{\rm err}/\sigma_{\kappa}$", fontsize=20)
if save_figs:
	plt.savefig("fig/ring_error.png")


plt.figure()
imgplot = plt.imshow(kappa_plane_orth_north.real,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("kappa plane orthographic projection north real")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_orth_north.png")

plt.figure()
imgplot = plt.imshow(kappa_plane_orth_north.real-kappa_orig_orth_north,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("kappa plane orthographic projection north error real")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_orth_north_error_real.png")

plt.figure()
imgplot = plt.imshow(kappa_plane_orth_north.imag,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("kappa plane orthographic projection north error imag")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_orth_north_error_imag.png")

plt.figure()
imgplot = plt.imshow(kappa_plane_orth_south.real,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("kappa plane orthographic projection south real")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_orth_south.png")

plt.figure()
imgplot = plt.imshow(kappa_plane_orth_south.real-kappa_orig_orth_south,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("kappa plane orthographic projection south error real")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_orth_south_error_real.png")

plt.figure()
imgplot = plt.imshow(kappa_plane_orth_south.imag,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("kappa plane orthographic projection south error imag")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_orth_south_error_imag.png")


# plt.figure()
# imgplot = plt.imshow(gamma_plane_orth_north.real,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.title("gamma plane orthographic projection north real")
# plt.axis('off')

# plt.figure()
# imgplot = plt.imshow(gamma_plane_orth_north.imag,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.title("gamma plane orthographic projection north imag")
# plt.axis('off')

# plt.figure()
# imgplot = plt.imshow(gamma_plane_orth_south.real,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.title("gamma plane orthographic projection south real")
# plt.axis('off')

# plt.figure()
# imgplot = plt.imshow(gamma_plane_orth_south.imag,interpolation='nearest')
# plt.colorbar(imgplot)
# plt.imshow(mask_north_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.title("gamma plane orthographic projection south imag")
# plt.axis('off')


print "start"
orth_proj_north, mask_north, orth_proj_south, mask_south = ssht.orthographic_projection(k_mw, L, resolution=500, Method=Method, rot=[0.0,np.pi/2,0.0])#zoom_region=np.pi/4)
print "finish"
plt.figure()
imgplot = plt.imshow(orth_proj_north,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_north, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("orthographic projection north")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_orig_orth_north.png")

plt.figure()
imgplot = plt.imshow(orth_proj_south,interpolation='nearest')
plt.colorbar(imgplot)
plt.imshow(mask_south, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.title("orthographic projection south")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_orig_orth_south.png")



print "start"
f_plot, masked_array = ssht.mollweide_projection(k_mw, L, resolution=500, Method=Method)
print "finish"
plt.figure()
imgplot = plt.imshow(f_plot,interpolation='nearest')
plt.colorbar(imgplot,fraction=0.025, pad=0.04)
plt.imshow(masked_array, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.gca().set_aspect("equal")
plt.title("Sphere kappa")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_orig_mollweide.png")


print "start"
f_real_plot, mask_real, f_imag_plot, mask_imag = ssht.mollweide_projection((k_mw-k_rec_mw), L, resolution=500, Method=Method)
print "finish"
plt.figure()
imgplot = plt.imshow(f_real_plot,interpolation='nearest')
plt.colorbar(imgplot,fraction=0.025, pad=0.04)
plt.imshow(mask_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.gca().set_aspect("equal")
plt.title("Sphere error real")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_mollweide_error_real.png")

plt.figure()
imgplot = plt.imshow(f_imag_plot,interpolation='nearest')
plt.colorbar(imgplot,fraction=0.025, pad=0.04)
plt.imshow(mask_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
plt.gca().set_aspect("equal")
plt.title("Sphere error imaginary")
plt.axis('off')
if save_figs:
	plt.savefig("fig/kappa_rec_mollweide_error_imag.png")

# print "start"
# f_plot, masked_array, dummy, dummy1 = ssht.mollweide_projection(kappa_plane, L, resolution=500, Method=Method)
# print "finish"
# plt.figure()
# imgplot = plt.imshow(f_plot,interpolation='nearest')
# plt.colorbar(imgplot,fraction=0.025, pad=0.04)
# plt.imshow(masked_array, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.gca().set_aspect("equal")
# plt.title("Plane kappa")
# plt.axis('off')

# print "start"
# f_real_plot, mask_real, f_imag_plot, mask_imag = ssht.mollweide_projection((k_mw-kappa_plane)*mask_small, L, resolution=500, Method=Method)
# print "finish"
# plt.figure()
# imgplot = plt.imshow(f_real_plot,interpolation='nearest')
# plt.colorbar(imgplot,fraction=0.025, pad=0.04)
# plt.imshow(mask_real, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.gca().set_aspect("equal")
# plt.axis('off')
# plt.title("Plane error real")
# plt.figure()
# imgplot = plt.imshow(f_imag_plot,interpolation='nearest')
# plt.colorbar(imgplot,fraction=0.025, pad=0.04)
# plt.imshow(mask_imag, interpolation='nearest', cmap=cm.gray, vmin=-1., vmax=1.)
# plt.gca().set_aspect("equal")
# plt.title("Plane error imaginary")
# plt.axis('off')


xticks = np.linspace(0.0,np.pi*2,5)
yticks = np.linspace(0.0,np.pi,5)

plt.figure()
imgplot = plt.imshow(kappa_plane.real,interpolation='nearest', extent=[0,2*np.pi,0,np.pi])
plt.colorbar(imgplot)
plt.title("Plane kappa")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\theta$")
if save_figs:
	plt.savefig("fig/kappa_rec_cylindrical.png")


plt.figure()
imgplot = plt.imshow((kappa_orig_plane-kappa_plane).real,interpolation='nearest', extent=[0,2*np.pi,0,np.pi])
plt.colorbar(imgplot)
plt.title("Plane error real")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\theta$")
if save_figs:
	plt.savefig("fig/kappa_rec_cylindrical_error_real.png")

plt.figure()
imgplot = plt.imshow((kappa_orig_plane-kappa_plane).imag,interpolation='nearest', extent=[0,2*np.pi,0,np.pi])
plt.colorbar(imgplot)
plt.title("Plane error imaginary")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\theta$")
if save_figs:
	plt.savefig("fig/kappa_rec_cylindrical_error_imag.png")

if show_figs:
	plt.show()

