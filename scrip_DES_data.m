clear all;
close all;


shear_real_data_construction

phi0 = L - 25 * (L/128);
theta0 = L/2;

figure(1)
subplot(121)
plot_orth_proj(real(g).*maskn,L, theta0, phi0,0)
title('real part of shear')
subplot(122)
plot_orth_proj(imag(g).*maskn,L,theta0, phi0,0)
title('imag part of shear')
sigma = pi /(3*180);

% Spherical KS convergence

k = shear2conv(g,L);
k = smooth_sph(k,sigma,L);

figure(2)
plot_orth_proj(k.*maskn,L,theta0, phi0,0)
title('spherical KS convergence')

% Plane KS convergence

k11 = shear2conv_plan(g,a,c,0,L);
k1 = zeros(L,2*L-1);
k1(a,c) = k11;
k1 = smooth_sph(k1,sigma,L);

figure(3)
plot_orth_proj(k1.*maskn,L,theta0, phi0,0)
title('plane KS convergence')

% Plane KS convergence with orthogonal projection

figure(4)
plot_orth_proj(g.*maskn,L,theta0, phi0,1)
title('plane KS convergence with orthognal projection')

% Plane KS convergence with plane DES data construction

figure(5)
plot_orth_proj_des
title('plane KS convergence with plane DES data construction')