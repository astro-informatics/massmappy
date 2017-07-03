import healpy as hp
import numpy as np
import pyssht as ssht
import cy_mass_mapping as mm
import matplotlib.pyplot as plt

def generate_k_maps(Cls,L,issimulation=False,ismasked=False,mask=None, Method='MW',sigma=None):
    if Method == 'MW':
        print "mm.generate_kappa_lm_mw"
        k_lm = mm.generate_kappa_lm_mw(np.array(Cls[:,1]), L)
        print "ssht.guassian_smoothing"
        ks_lm = ssht.guassian_smoothing(k_lm, L, sigma_in=np.pi/256)
        print "ssht.inverse"
        k = ssht.inverse(ks_lm, L, Reality=True, Method=Method)
        print "mm.kappa_lm_to_gamma_lm_mw"
        gamma_lm = mm.kappa_lm_to_gamma_lm_mw(ks_lm, L)
        print "ssht.inverse"
        gamma = ssht.inverse(gamma_lm, L, Method=Method, Spin=2)
        #add noise AFTER shear is calculated, not before
        if issimulation == True:
            shear = gamma/(1.0-k) + (gamma.std()*(np.random.randn(L,2*L-1) + 1j*np.random.randn(L,2*L-1)))
            #shear = gamma/(1.0-k) +  (0.01*(np.random.randn(L,2*L-1) + 1j*np.random.randn(L,2*L-1)))
        else:
            shear = gamma/(1.0-k)
        print "Reconstruct on the sphere MW"
        print "sigma= ",sigma
        if sigma is not None:
            if Method == 'MW':
                if ismasked == True:
                    k_rec = mm.reduced_shear_to_kappa_mw(shear*mask, L, Method=Method, sigma=sigma,tol_error=1E-1)
                else:
                    k_rec = mm.reduced_shear_to_kappa_mw(shear, L, Method=Method, sigma=sigma,tol_error=1E-1)
        else:
            if Method == 'MW':
                if ismasked == True:
                    k_rec = mm.reduced_shear_to_kappa_mw(shear*mask, L, Method=Method, tol_error=1E-1)
                else:
                    k_rec = mm.reduced_shear_to_kappa_mw(shear, L, Method=Method, tol_error=1E-1)
    elif Method == 'HP':
        print "mm.generate_kappa_lm_hp"
        k_lm = mm.generate_kappa_lm_hp(np.array(Cls[:,1]), L)
        print "hp.sphtfunc.smoothalm"
        ks_lm = hp.sphtfunc.smoothalm(k_lm,sigma=np.pi/256)
        print "hp.alm2map"
        #want nside to give similar resolution of equialent L, make sure L is a power of 2, ie 2^n. why does alm2map require lmax,mmax=L-1?
        k = hp.alm2map(ks_lm,nside=L,lmax=L-1,mmax=L-1)
        print "mm.kappa_lm_to_gamma_lm_hp"
        gamma_lm = mm.kappa_lm_to_gamma_lm_hp(ks_lm, L)[0]
        print "hp.alm2map"
        gamma = hp.alm2map(gamma_lm, nside=L)
        #add noise AFTER shear is calculated, not before
        if issimulation == True:
            shear = gamma/(1.0-k) + 0.01*(np.random.randn(gamma.shape[0]) + 1j*np.random.randn(gamma.shape[0]))
        else:
            shear = gamma/(1.0-k)
        if ismasked == True:
            #NOTE: currently cannot do masks in HP
            [k_rec_E,k_rec_B] = mmhp.reduced_shear_to_kappa_hp(np.ascontiguousarray(np.real(shear)*mask),np.ascontiguousarray(np.imag(shear)*mask), L,L, tol_error=1E-1)
            k_rec = k_rec_E + k_rec_B*1j
        else:
            [k_rec_E,k_rec_B] = mmhp.reduced_shear_to_kappa_hp(np.ascontiguousarray(np.real(shear)),np.ascontiguousarray(np.imag(shear)), L,L, tol_error=1E-1)
            k_rec = k_rec_E + k_rec_B*1j
    return [k_lm, ks_lm,k, gamma_lm, gamma, shear, k_rec]

def gradient_weights(L, Method="MW"):

    thetas, phis = ssht.sample_positions(L, Method=Method, Grid=True)

    theta_weights = np.ones(thetas.shape)*L/np.pi
    phi_weights = np.ones(thetas.shape)*L/np.pi
    phi_weights[0:thetas.shape[0]-1,:]   = L/(np.pi*np.sin(thetas[0:thetas.shape[0]-1,:]))

    return theta_weights, phi_weights

def take_gradient(f, L, Method="MW"):

    n_theta, n_phi = ssht.sample_shape(L, Method=Method)
    theta_weights, phi_weights = gradient_weights(L, Method=Method)

    g_theta = np.zeros((n_theta,n_phi))
    g_phi   = np.zeros((n_theta,n_phi))

    g_theta[:n_theta-1,:] = f[1:,:] - f[:n_theta-1,:]
    g_theta[n_theta-1,:]  = 0

    g_phi[:,:n_phi-1]     = f[:,1:]-f[:,:n_phi-1]
    g_phi[:,n_phi-1]=0

    return theta_weights*g_theta, g_phi*phi_weights

L = 256
Method = "MW"
sigma_scale = 5.0
sigma_MW=2.*20.0*np.pi/(60.*180.*2.355)*sigma_scale
Cls = np.loadtxt("data/cls_ap.txt")
issimulation = True
looparound = False
nu_step = 0.25
isGTzero = True
ismasked = False
mask = np.zeros((L,2*L-1), dtype=float)
Reality = True
Spin = 0

[k_lm, ks_lm, k, gamma_lm, gamma, shear, k_rec] = generate_k_maps(Cls,L,issimulation, ismasked, mask, Method, sigma=sigma_MW)

nu = (np.real(k_rec)-np.nanmean(np.real(k_rec)))/np.real(np.nanstd(k_rec))

amap = nu

nside_match = int(np.ceil(np.sqrt(L*(2*L-1)/12)))
nside = 2**int(np.ceil(np.log2(nside_match)))

plt.imshow(amap)

alm = ssht.forward(amap, L, Reality=Reality, Method=Method, Spin=Spin)
alm_hp = mm.lm2lm_hp(alm,L)
amap_hp, d_theta, d_phi = hp.sphtfunc.alm2map_der1(alm_hp,nside)
# hp.visufunc.mollview(amap_hp)

d_theta_lm = hp.sphtfunc.map2alm(d_theta,lmax=L-1)
d_phi_lm = hp.sphtfunc.map2alm(d_phi,lmax=L-1)
d_theta_mw_lm = mm.lm_hp2lm(d_theta_lm.astype('complex'),L)
d_phi_mw_lm = mm.lm_hp2lm(d_phi_lm.astype('complex'),L)

d_theta_mw = ssht.inverse(d_theta_mw_lm,L,Spin=Spin)
d_phi_mw = ssht.inverse(d_phi_mw_lm,L,Spin=Spin)

plt.figure()
plt.imshow(d_theta_mw.real)
plt.colorbar()
plt.figure()
plt.imshow(d_phi_mw.real)
plt.colorbar()

g_theta, g_phi = take_gradient(amap,L)
plt.figure()
plt.imshow(g_theta)
plt.colorbar()
plt.figure()
plt.imshow(g_phi)
plt.colorbar()

print d_theta_mw.real/g_theta
plt.show()

