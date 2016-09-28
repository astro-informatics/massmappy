function [ out_f ] = smooth_sph( f,sigma,L )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out_f = ssht_inverse(smooth_lm(ssht_forward(f,L),sigma,L),L);

end

