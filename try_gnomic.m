
L = 256;
addpath('/unsafe1/opt/ssht/src/matlab')
% 
% m = rotation_mask(mask,1,0,0,256);
% f = zeros(256,511);
% f(40:50,40:50) = ones(11,11);
% 
% 
% [a,b,c] = angle2notrhpole(45,45,L);
% 
% figure(1)
% subplot(121)
% ssht_plot_mollweide(f,256);
% subplot(122)
% ssht_plot_mollweide(rotation_sphere(f,a,b,c,256),256);



% MASK WITH INTERESSANT REGION ON THE EQUATOR :
% rotation_mask(mask,0,2*pi/7-0.2,0.3 + pi/2,256) -> rotation of angle
% -pi/2 to obtain that at north pole.

% figure(3)
% subplot(221)
% ssht_plot_mollweide(rotation_sphere(g,0,pi,pi+0.5,64),64);
% subplot(222)
% ssht_plot_mollweide(rotation_sphere(g,0,0,0,64),64);
% subplot(223)
% ssht_plot_mollweide(rotation_sphere(g,0,0.2,0,64),64);
% subplot(224)
% ssht_plot_mollweide(rotation_sphere(g,0,pi,0,64),64);

[a,b,c] = angle2notrhpole(m1,m2,256);
figure(4) 
subplot(111)
plot_gnomonic_proj(g,256,a,b,c,mask);
title('Shear')
% rotation to center it : plot_gnomonic_proj(mask,256,0.4,pi/2,0.4,mask) 
% figure(4)
% subplot(221)
% plot_gnomonic_proj(rotation_sphere(g,0,pi/2,pi/2,256),256);
% subplot(223)
% plot_gnomonic_proj(rotation_sphere(g,0,pi/2,pi/2+0.3,256),256);
% subplot(224)
% plot_gnomonic_proj(rotation_sphere(g,0,pi/2,pi/2+0.5,256),256);
% subplot(224)
% plot_gnomonic_proj(rotation_sphere(g,0.5,pi/2,pi/2+0.15,256),256);