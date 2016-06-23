function [ k_l ] = low_res_conv( k_h, L )

addpath('/unsafe1/opt/ssht/src/matlab')

% k_h resolution is 4*L, and k_l resolution is L

klm = ssht_forward(k_h,4*L, 'Reality', true);

klm = smooth_lm(klm, 0.01, 4*L);

k_l = ssht_inverse(klm(1:L^2),L, 'Reality', true);

end

