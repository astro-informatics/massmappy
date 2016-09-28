
addpath('/unsafe1/opt/ssht/src/matlab')


rng(3)

N = 100;
L = 2048;
T1 = L * (2*L - 1);

D = zeros(L,1);

klm = [];
gammalm= [];

for l = 1:L
    D(l+1) = -1 /(l*(l+1)) * sqrt((l+2)*(l+1)*l*(l-1));
end

for l = 0:L-1  
    
    km = zeros(2*l+1,1);
    km(l+1:2*l+1) = normrnd(0, sqrt( Cell(l+1) ), [ l+1 ,1]) + 1i * normrnd(0, sqrt( Cell(l+1) ), [ l+1 ,1]) ;
    
    for m = 1:l
       km(l-m+1) = (-1)^m * km(l+m);
    end
    
    klm = [klm; km ];
    
    gammalm = [gammalm; D(l+1)*km];
    
end

klm(1:4) = [0;0;0;0];

k = real( ssht_inverse(klm,L) );

gamma = ssht_inverse(gammalm,L,'Spin',2);

M = mask;

gamma = gamma .* M ; 

[k1_1,k2_1] = shear2conv_plan(gamma,a,c,0,L);
k1 = zeros(L,2*L-1);
k1(a,c) = k1_1;
k1s = smooth_sph(k1,pi/(180*3),L);
ks = smooth_sph(k,pi/(180*3),L);

k2 = zeros(L,2*L-1);
k2(a,c) = k1_1;
k2s = smooth_sph(k2,pi/(180*3),L);
ks = smooth_sph(k,pi/(180*3),L);

figure(1)
subplot(121)
plot_orth_proj(ks.*maskn,L,L/2,L-400,0);
title('Simulated convergence');
caxis([min(min(real(ks))), max(max(real(ks)))])
subplot(122)
plot_orth_proj(k1s.*maskn,L,L/2,L-400,0);
title('plane convergence k1 (sans normalisation)');
caxis([min(min(real(ks))), max(max(real(ks)))])

figure(2)
subplot(111)
plot_orth_proj(abs(real(ks)-real(k1s)).*maskn,L,L/2,L-400,0);
% caxis([min(min(abs(real(ks)-real(k1s)))), max(max(abs(real(ks)-real(k1s))))])
title('Simulated convergence');



