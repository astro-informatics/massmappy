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
Cell = ones(1,2048);

for l = 0:L-1  
    
    km = zeros(2*l+1,1);
    km = normrnd(0, sqrt( Cell(l+1) ), [ l ,1]) + 1i * normrnd(0, sqrt( Cell(l+1) ), [ l ,1]) ;
    km(l+1,1) = normrnd(0, sqrt( Cell(l+1) ));
    
    for m = 1:l
        km(l+1+m,1) = (-1).^(m) * conj(km(l+1-m,1)); 
    end
    klm = [klm; km ];
    
    gammalm = [gammalm; D(l+1)*km];
    
end

klm(1:4) = [0;0;0;0];

k = ssht_inverse(klm,L);

gamma = ssht_inverse(gammalm,L,'Spin',2);


% m = rotation_mask(mask,0,2*pi/7+0.2,0,256);
% noise = realistic_noise(mask, L,2*L-1, L);

N = 50;
N1 = L/2-N/2+1  : L/2+N/2 ;
N2 = L-N/2+1 : L+N/2;

mask = zeros(L,2*L-1);
mask(N1, N2) = ones(size(mask(N1, N2)));
gamma1 = gamma.*mask;

k1 = shear2conv(gamma1,L);

% 
% N = 100;
% th = (L-N)/2 +1 : (L+N)/2;
% ph = L - N/2 + 1: L + N/2;
% maskn = mask2nan(mask);

N = 50;
N1 = L/2-N/2+1  : L/2+N/2 ;
N2 = L-N/2+1 : L+N/2;

figure(1)
subplot(222)
imagesc(real(k(N1,N2)));
title('Simulated Convergence')

subplot(224)
imagesc(real(k1(N1,N2)));
title('Convergence k1 from noisy shear')

subplot(221)
imagesc(real(gamma(N1,N2)));
title('Corresponding shear')

subplot(223)
imagesc(abs(k(N1,N2)-real(k1(N1,N2))));
title('diff k k1')

% figure(2)
% subplot(211)
% plot_orth_proj(abs(k-k1).*maskn,L,L/2,L-100);
% title('abs(k-k1)')
% subplot(212)
% plot_orth_proj(abs(gamma-gamma1).*maskn,L,L/2,L-100);
% title('abs(g-g1)')
% % 
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