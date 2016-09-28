function [ out_m ] = smooth_mask( m , L)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
addpath('/unsafe1/opt/ssht/src/matlab')

m_lm = ssht_forward(m,L,'Reality', true);
m_lm = smooth_lm(m_lm,0.01,L);
out_m = ssht_inverse(m_lm,L,'Reality', true);

out_m = double( out_m > 0.3 );

end

