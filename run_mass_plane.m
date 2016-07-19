close all; 
clear all;

read_cell_data;

% L = 200:50:500;
% L = [L,600 :200: 2000];
% Cell = ones(1,2200);


% 
% theta0 = 0:50:400;
% N = 50:50:450;
% 
r = 0.25 : 0.25 : 3;
T = size(r,2);
N= 80;
L=1024;
error_p = zeros(1,T);
error_s = zeros(1,T);
rms_k = zeros(1,T);

for l = 1:T;
     
    l
    [s,p,k] = mass_map_plane_matrix_function(Cell, L,N,2,0,r(l));
    error_p(l) = p;
    error_s(l) = s;
    rms_k(l) = k;
     
end

 ang =  N./ L * 180;
% dec = theta0 * 180/L;
figure(3)
plot(r,error_p./rms_k, 'o-',r,error_s./rms_k,'o-');
legend('error plane KS', 'error spherical KS')
title('error of KS, no smoothing, bufferRatio = 2')
xlabel('rectangle ration')
ylabel('error percentage (RMS)')
% figure(1)
% plot(r,rms_k, 'o-');
% legend('error plane KS', 'error spherical KS')
% title('rms of simulated convergence k, no smoothing, bufferRatio = 1')
% xlabel('rect_ratio')
% ylabel('RMS k')
% 
% figure(3)
% subplot(111)
% plot(rms_k,error_p, 'o');
% legend('error plane KS', 'error spherical KS')
% title('rms of simulated convergence k, no smoothing, bufferRatio = 1')
% xlabel('RMS k')
% ylabel('RMS plane error')