function [ out_n ] = realistic_noise(mask, t1, t2, L )

out_n = zeros(t1, t2);

Ngal = 10^9;
M = sum(sum(mask));
dThe = pi/(L-1);

k = 1;
l = 1;
dPh2 = 0;
Norm = 6.73e-05;
while l <= t1
    
    k=1;
    dPh1 = dPh2;
    dPh2 = sin(pi*l/L)*2*pi/(2*L-1);
    dOm = abs( dThe*(dPh2+dPh1)/2 )/Norm;
    ngal = Ngal/M*dOm

    while k <= t2
       
        if(mask(l,k) ~= 0)

           out_n(l,k) = rayleigh_distrib(floor(ngal));
           
           k = k+1;
        else
          k = k+1;
        end
        k;
    end
    l = l+1
end
end

