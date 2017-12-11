import numpy as np
import matplotlib.pyplot as plt
import pyssht as ssht
import massmappy.cy_mass_mapping as mm
from matplotlib import cm

show_figs = False
save_figs = True

Cls = np.loadtxt("data/cls_ap.txt")

# error_cylindrical_E = np.loadtxt("data/high_res_results/keep_2/error_cylindrical_E.txt")
# error_mercator_E = np.loadtxt("data/high_res_results/keep_2/error_mercater_E.txt")
# error_sine_E = np.loadtxt("data/high_res_results/keep_2/error_sine_E.txt")
# error_orthographic_E = np.loadtxt("data/high_res_results/keep_2/error_orthographic_E.txt")
# error_steriographic_E = np.loadtxt("data/high_res_results/keep_2/error_steriographic_E.txt")
# error_gnomic_E = np.loadtxt("data/high_res_results/keep_2/error_gnomic_E.txt")
# error_cylindrical_B = np.loadtxt("data/high_res_results/keep_2/error_cylindrical_B.txt")
# error_mercator_B = np.loadtxt("data/high_res_results/keep_2/error_mercater_B.txt")
# error_sine_B = np.loadtxt("data/high_res_results/keep_2/error_sine_B.txt")
# error_orthographic_B = np.loadtxt("data/high_res_results/keep_2/error_orthographic_B.txt")
# error_steriographic_B = np.loadtxt("data/high_res_results/keep_2/error_steriographic_B.txt")
# error_gnomic_B = np.loadtxt("data/high_res_results/keep_2/error_gnomic_B.txt")

error_cylindrical_E = np.loadtxt("data/high_res_results/error_cylindrical_E.txt")
error_mercator_E = np.loadtxt("data/high_res_results/error_mercater_E.txt")
error_sine_E = np.loadtxt("data/high_res_results/error_sine_E.txt")
error_orthographic_E = np.loadtxt("data/high_res_results/error_orthographic_E.txt")
error_steriographic_E = np.loadtxt("data/high_res_results/error_steriographic_E.txt")
error_gnomic_E = np.loadtxt("data/high_res_results/error_gnomic_E.txt")
error_cylindrical_B = np.loadtxt("data/high_res_results/error_cylindrical_B.txt")
error_mercator_B = np.loadtxt("data/high_res_results/error_mercater_B.txt")
error_sine_B = np.loadtxt("data/high_res_results/error_sine_B.txt")
error_orthographic_B = np.loadtxt("data/high_res_results/error_orthographic_B.txt")
error_steriographic_B = np.loadtxt("data/high_res_results/error_steriographic_B.txt")
error_gnomic_B = np.loadtxt("data/high_res_results/error_gnomic_B.txt")


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

N_angle = 12*2000/500
angle, angle_step = np.linspace(0.0, np.pi/2, N_angle, retstep=True)

CFHTLens = np.sqrt(154.)
DES_SV = np.sqrt(139.)
KiDS = np.sqrt(1500.)
DES_full = np.sqrt(5000.)
Euclid = np.sqrt(15000.)
LSST = np.sqrt(15000.)

label_buffer = 2.0
ell = 200
#normalisation = np.sqrt(Cls[ell,1]*ell*(ell+1)/(2*np.pi))
normalisation = 1.0


plt.figure()
plt.plot([DES_SV,DES_SV], [0,1], 'k--')
plt.text(DES_SV+label_buffer, 0.5, "CFHTLens\nDES SV", fontsize=18)
plt.plot([DES_full,DES_full], [0,1], 'k--')
plt.text(DES_full+label_buffer, 0.9, "DES full", fontsize=18)
plt.plot([LSST, LSST], [0,1], 'k--')
plt.text(LSST+label_buffer, 0.85, "Euclid\nLSST", fontsize=18)
plt.errorbar(angle*2*180/np.pi, error_cylindrical_E_av/normalisation,color='#377eb8',\
	yerr=error_cylindrical_E_std/normalisation, label="Cylindrical")
plt.errorbar(angle*2*180/np.pi, error_mercator_E_av/normalisation,color='#ff7f00',\
	yerr=error_mercator_E_std/normalisation, label="Mercator")
plt.errorbar(angle*2*180/np.pi, error_sine_E_av/normalisation,color='#4daf4a',\
	yerr=error_sine_E_std/normalisation, label="Sinusiodal")
plt.errorbar(angle*2*180/np.pi, error_orthographic_E_av/normalisation,color='#f781bf',\
	yerr=error_orthographic_E_std/normalisation, label="Orthographic")
plt.errorbar(angle*2*180/np.pi, error_steriographic_E_av/normalisation,color='#a65628',\
	yerr=error_steriographic_E_std/normalisation, label="Stereographic")
plt.errorbar(angle*2*180/np.pi, error_gnomic_E_av/normalisation,color='#984ea3',\
	yerr=error_gnomic_E_std/normalisation, label="Gnomonic")
plt.xlabel(r"$\Theta \, ({\rm degrees})$", fontsize=20)
plt.ylabel(r"$\sigma_{\rm err}/\sigma_{\kappa}$", fontsize=20)
plt.legend(loc=2)
plt.xlim([0,180])
plt.ylim([0,1.0])
if save_figs:
	plt.savefig("fig/ring_error_E.pdf")

plt.figure()
plt.plot([DES_SV,DES_SV], [0,1], 'k--')
plt.text(DES_SV+label_buffer, 0.5, "CFHTLens\nDES SV", fontsize=18)
plt.plot([DES_full,DES_full], [0,1], 'k--')
plt.text(DES_full+label_buffer, 0.9, "DES full", fontsize=18)
plt.plot([LSST, LSST], [0,1], 'k--')
plt.text(LSST+label_buffer, 0.85, "Euclid\nLSST", fontsize=18)
plt.errorbar(angle*2*180/np.pi, error_cylindrical_B_av/normalisation,color='#377eb8',\
	yerr=error_cylindrical_B_std/normalisation, label="Cylindrical")
plt.errorbar(angle*2*180/np.pi, error_mercator_B_av/normalisation,color='#ff7f00',\
	yerr=error_mercator_B_std/normalisation, label="Mercator")
plt.errorbar(angle*2*180/np.pi, error_sine_B_av/normalisation,color='#4daf4a',\
	yerr=error_sine_B_std/normalisation, label="Sinusiodal")
plt.errorbar(angle*2*180/np.pi, error_orthographic_B_av/normalisation,color='#f781bf',\
	yerr=error_orthographic_B_std/normalisation, label="Orthographic")
plt.errorbar(angle*2*180/np.pi, error_steriographic_B_av/normalisation,color='#a65628',\
	yerr=error_steriographic_B_std/normalisation, label="Stereographic")
plt.errorbar(angle*2*180/np.pi, error_gnomic_B_av/normalisation,color='#984ea3',\
	yerr=error_gnomic_B_std/normalisation, label="Gnomonic")
plt.xlabel(r"$\Theta \, ({\rm degrees})$", fontsize=20)
plt.ylabel(r"$\sigma_{\rm err}/\sigma_{\kappa}$", fontsize=20)
plt.legend(loc=2)
plt.xlim([0,180])
plt.ylim([0,1.0])
if save_figs:
	plt.savefig("fig/ring_error_B.pdf")

if show_figs:
	plt.show()
