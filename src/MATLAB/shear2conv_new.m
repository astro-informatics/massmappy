function [ out_k] = shear2conv_new( gamma, L )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%g = gammalm;

gamma_lm = ssht_forward(gamma,L,'Spin',2);

for el = 2:L-1
    for m = -el:el
    
        ind = ssht_elm2ind(el, m);      
        kappa_lm(ind) = - gamma_lm(ind) .* el .* (el+1) ./ ...
            sqrt((el+2).*(el+1).*el.*(el-1));
        
    end
end

out_k = ssht_inverse(kappa_lm, L); 
%out_k = ssht_inverse(g, L,'Spin',2); 
end

