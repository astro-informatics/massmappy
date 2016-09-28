close all; 
clear all;

read_cell_data;

test = 1; % 1: L, 2: rectangle ratio, 3: buffer ratio
band_limit_mask = 1; % 0: not bandlimited, 1: band limited, 2: appodised mask
mask_type = 1; % 0: square/rectangle mask, 1: cirlular mask
N= 100;
L=1024;

if test == 1,
 L = 100:50:500;
T = size(L,2);
end

% 
% theta0 = 0:50:400;
% N = 50:50:450;
% 
if test == 2,
r = 0.25 : 0.25 : 3.0;
T = size(r,2);
end
%
if test == 3,
br = 1.0 : 1.0 : 5.0;
T = size(br,2);
end

error_p = zeros(1,T);
error_s = zeros(1,T);
rms_k = zeros(1,T);

for l = 1:T;
    l
    if test == 1,
        if 2.0*N > L(l),
            br = L(l)/N
        else
            br = 2.0
        end
        [s,p,k] = mass_map_plane_matrix_function(Cell, L(l),N,br,0,1.0, band_limit_mask, mask_type);
    end
    if test == 2,
        [s,p,k] = mass_map_plane_matrix_function(Cell, L,N,2.0,0,r(l), band_limit_mask, mask_type);
    end
    if test == 3,
        [s,p,k] = mass_map_plane_matrix_function(Cell, L,N,br(l),0,1.0, band_limit_mask, mask_type);
    end
    error_p(l) = p;
    error_s(l) = s;
    rms_k(l) = k;
     
end

 ang =  N./ L * 180;
% dec = theta0 * 180/L;
if test == 1,
figure(2)
plot(ang,error_p./rms_k, 'o-',ang,error_s./rms_k,'o-');
legend('error plane KS', 'error spherical KS')
title_string = 'error of KS, no smoothing, bufferRatio = 2'
if band_limit_mask == 0,
title_string = strcat(title_string,', not band limited mask');
else
title_string = strcat(title_string, ', band limited mask');
end
if mask_type == 0,
title_string = strcat(title_string,', square mask');
else
title_string = strcat(title_string, ', circular mask');
end
title(title_string);
xlabel('angle (degrees)')
ylabel('error percentage (RMS)')
end

if test == 2,
figure(2)
plot(r,error_p./rms_k, 'o-',r,error_s./rms_k,'o-');
legend('error plane KS', 'error spherical KS')
title('error of KS, no smoothing, bufferRatio = 2')
xlabel('rectangle ratio')
ylabel('error percentage (RMS)')
end
 %  
 
 
if test == 3,
     figure(2)
plot(br,error_p./rms_k, 'o-',br,error_s./rms_k,'o-');
legend('error plane KS', 'error spherical KS')
title('error of KS, no smoothing')
xlabel('buffer ratio')
ylabel('error percentage (RMS)')
end
%  figure(1)
%  plot(br,rms_k, 'o-');
%  legend('error plane KS', 'error spherical KS')
%  title('rms of simulated convergence k, no smoothing')
%  xlabel('buffer ratio')
%  ylabel('RMS k')
% % 
%  figure(3)
%  subplot(111)
%  plot(rms_k,error_p, 'o');
%  legend('error plane KS', 'error spherical KS')
%  title('rms of simulated convergence k, no smoothing')
%  xlabel('RMS k')
%  ylabel('RMS plane error')