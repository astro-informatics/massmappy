import numpy as np
import pyssht as ssht

def gradient_weights(L, Method="MW"):

	thetas, phis = ssht.sample_positions(L, Method=Method, Grid=True)

	theta_weights = np.ones(thetas.shape)
	phi_weights   = 1.0/np.sin(thetas)

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

	return theta_weights*g_theta, phi_weights*g_phi

L = 4

print gradient_weights(L)

n_theta, n_phi = ssht.sample_shape(L)
f = np.ones((n_theta,n_phi))

print take_gradient(f, L)


