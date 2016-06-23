function [ out_k ] = reducedshear2conv( g, L )

addpath('/unsafe1/opt/ssht/src/matlab')
eps = 1e-04;

out_k = shear2conv(g, L);
k1 = zeros(size(g));
gamma = g;

n = 1;

while( sum( sum( abs(k1-out_k) ) ) > eps )
    
    gamma = g .* (1-out_k);
    
    k1 = out_k;
    
    out_k = shear2conv(gamma,L);
    
    n = n + 1
    
end


end

