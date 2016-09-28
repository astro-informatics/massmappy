function [ out_flm ] = smooth_lm( flm, sigma, L )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

out_flm = zeros(size(flm));
smoo = exp(-( sigma*(0:L-1) ).^2);

for k = 1:size(flm,1)
    
    [l,m] = ssht_ind2elm(k);
    out_flm(k) = smoo(l+1)*flm(k);
    
end

end

