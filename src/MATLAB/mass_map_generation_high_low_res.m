
addpath('/unsafe1/opt/ssht/src/matlab')


rng(11)

L = 512;
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

k = ssht_inverse(klm,L, 'Reality', true);

gamma = ssht_inverse(gammalm,L,'Spin',2);

L = 128;

%noise = realistic_noise(mask, L,2*L-1, L);
noise = zeros(L,2*L-1);

% gamma1 = (gamma + noise).*mask;
gamma1 = low_shear_res(gamma,L);
k = low_res_conv(k,L);

% k = low_shear_res(k,L);
% k = ssht_inverse( smooth_lm(ssht_forward(k1,L,'Reality',true),0.01,L),L,'Reality', true);

k1 = shear2conv(gamma1,L);
k1 = ssht_inverse( smooth_lm(ssht_forward(k1,L,'Reality',true),0.01,L),L,'Reality', true);


figure(1)
subplot(222)
ssht_plot_mollweide(k,L);
title('Simulated Convergence')

subplot(224)
ssht_plot_mollweide(k1,L);
title('Convergence k1 from low resolution shear')

subplot(221)
ssht_plot_mollweide(abs(gamma),4*L);
title('Corresponding shear')

subplot(223)
ssht_plot_mollweide(abs(gamma1),L);
title('low resolution shear')

figure(2)
subplot(111)
ssht_plot_mollweide(abs(k-k1),L);
title('abs(k-k1)')

