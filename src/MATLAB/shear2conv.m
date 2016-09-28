function [ out_k] = shear2conv( gamma, L )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

addpath('/unsafe1/opt/ssht/src/matlab')

D = zeros(L,1);

for l = 1:L
    D(l+1) = -1 /(l*(l+1)) * sqrt((l+2)*(l+1)*l*(l-1));  % l'indice dans D qui correspond à l est égal à (l+1) car l commence à 0.
end

%g = gammalm;

g = ssht_forward(gamma,L,'Spin',2);

d = zeros(size(g));

d(1) = 1;

for ind = 2:size(g,1)  
    
    [l,m] = ssht_ind2elm(ind);

    d(ind) = D(l+1);
    if(d(ind) == 0)   % d(ind) = 0 signifie l = 1 et donc gammalm(ind) = 0, donc klm(ind) = 0.
        d(ind) = 1;
    end
end

klm = g./d;

out_k = ssht_inverse(klm, L); 
%out_k = ssht_inverse(g, L,'Spin',2); 
end

