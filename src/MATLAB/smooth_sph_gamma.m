function [ out_f ] = smooth_sph_gamma( f,sigma,L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out_f = ssht_inverse(smooth_lm(ssht_forward(f,L,'Spin',2),sigma,L),L,'Spin',2);

end

