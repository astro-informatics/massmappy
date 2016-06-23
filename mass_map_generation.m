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
    
    km = normrnd(0, sqrt(Cell(l+1) ), [ 2*l+1 ,1]) + 1i * normrnd(0, sqrt(Cell(l+1) ), [ 2*l+1 ,1]) ;
  
    klm = [klm; km ];
    
    gammalm = [gammalm; D(l+1)*km];
    
end

klm(1:4) = [0;0;0;0];

k = ssht_inverse(smooth_lm(klm,1/(3*180)*pi,L),L, 'Reality', true);

gamma = ssht_inverse(smooth_lm(gammalm,1/(3*180)*pi,L),L,'Spin',2);


% m = rotation_mask(mask,0,2*pi/7+0.2,0,256);
% noise = realistic_noise(mask, L,2*L-1, L);

gamma1 = gamma.*mask;

k1 = shear2conv(gamma1,L);


N = 100;
th = (L-N)/2 +1 : (L+N)/2;
ph = L - N/2 + 1: L + N/2;
maskn = mask2nan(mask);


figure(1)
subplot(222)
plot_orth_proj(k.*maskn,L,L/2,L-100);
title('Simulated Convergence')

subplot(224)
plot_orth_proj(k1.*maskn,L,L/2,L-100);
title('Convergence k1 from noisy shear')

subplot(221)
plot_orth_proj(gamma.*maskn,L,L/2,L-100);
title('Corresponding shear')

subplot(223)
plot_orth_proj(gamma1.*maskn,L,L/2,L-100);
title('Noisy shear')

figure(2)
subplot(211)
plot_orth_proj(abs(k-k1).*maskn,L,L/2,L-100);
title('abs(k-k1)')
subplot(212)
plot_orth_proj(abs(gamma-gamma1).*maskn,L,L/2,L-100);
title('abs(g-g1)')
% 
% figure(3)
% subplot(222)
% b = imagesc(k(th, ph).*maskn(th,ph));
% set(b,'AlphaData',~isnan(maskn(th,ph)));
% title('Simulated Convergence')
% 
% subplot(224)
% b = imagesc(k1(th, ph).*maskn(th,ph));
% set(b,'AlphaData',~isnan(maskn(th,ph)));
% title('Convergence k1 from noisy shear')
% 
% subplot(221)
% b = imagesc(abs(gamma(th, ph)).*maskn(th,ph));
% set(b,'AlphaData',~isnan(maskn(th,ph)));
% title('Corresponding shear')
% 
% subplot(223)
% b = imagesc(abs(gamma1(th, ph)).*maskn(th,ph));
% set(b,'AlphaData',~isnan(maskn(th,ph)));
% title('Noisy shear')
% 
% figure(4)
% subplot(211)
% imagesc(abs(k(th, ph)-k1(th, ph)));
% title('abs(k-k1)')
% subplot(212)
% imagesc(abs(gamma(th, ph)-gamma1(th, ph)));
% title('abs(g-g1)')