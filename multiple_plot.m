% k11 = ssht_inverse(smooth_lm(ssht_forward(k1,L, 'Reality', true),0.005,L),L, 'Reality', true);

figure(1)
subplot(111)
imagesc(real(ksp));
title('Simulated convergence k')

figure(2)
subplot(111)
imagesc(real(k1sp));
title('k1 obtained using spherical KS')

figure(3)
subplot(111)
imagesc(real(k2p));
title('k2 obtained using spherical KS')


figure(4)
subplot(111)
imagesc(abs(real(ksp)-real(k1sp)))
title('k-k1');

figure(5)
subplot(111)
imagesc(abs(real(k1sp)-real(k2p)))
title('k1-k2');
