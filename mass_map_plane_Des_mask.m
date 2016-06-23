
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
    
    km = normrnd(0, sqrt( Cell(l+1) ), [ 2*l+1 ,1]) + 1i * normrnd(0, sqrt( Cell(l+1) ), [ 2*l+1 ,1]) ;
    
    klm = [klm; km ];
    
    gammalm = [gammalm; D(l+1)*km];
    
end

klm(1:4) = [0;0;0;0];

k = real( ssht_inverse(klm,L) );

gamma = ssht_inverse(gammalm,L,'Spin',2);

M = mask;

gamma = gamma .* M ; 

% k1_1 = shear2conv_plan(gamma,a,c,12*pi/(180*3),L);
% k1 = zeros(L,2*L-1);
% k1(a,c) = k1_1;

figure(5)
subplot(121)
plot_orth_proj(smooth_sph(k,pi/(180*3),L).*maskn,L,L/2,L-400,0);
title('Simulated convergence');

subplot(122)
plot_orth_proj(gamma.*maskn,L,L/2,L-400,1);
title('plane convergence k1');

