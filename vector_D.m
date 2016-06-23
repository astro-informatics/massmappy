function [ out_D ] = vector_D(L)
%calcul de la matrice d qui contient les valeurs de D(l) associées à k(l,m)
%pour m entre -l et l;

T = ssht_elm2ind(L-1,L-1);

d = zeros(L,1);
for l = 1:L
    d(l+1) = -1 /(l*(l+1)) * sqrt((l+2)*(l+1)*l*(l-1));
end

out_D = zeros(T,1);

out_D(1) = 1;

for ind = 2:T
    
    [l,m] = ssht_ind2elm(ind);

    out_D(ind) = d(l+1);
    if(out_D(ind) == 0)   % d(ind) = 0 signifie l = 1 et donc gammalm(ind) = 0, donc klm(ind) = 0.
        out_D(ind) = 1;
    end
end


end

