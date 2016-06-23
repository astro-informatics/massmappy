k11 = ssht_inverse(smooth_lm(ssht_forward(k1,L, 'Reality', true),0.005,L),L, 'Reality', true);

figure(1)
subplot(111)
plot_orth_proj(real(gamma).*maskn,L,L/2,L-100);
title('Real part of shear g obtained from k')

figure(8)
subplot(111)
plot_orth_proj(real(gamma1).*maskn,L,L/2,L-100);
title('Observed shear g1, with Euclid mask and noise added')

figure(2)
subplot(111)
plot_orth_proj(k.*maskn,L,L/2,L-100);
title('Simulated Convergence k')


figure(9)

subplot(111)
plot_orth_proj(k1.*maskn,L,L/2,L-100);
title('Convergence k1 computed from g1')
colorbar

figure(3)
subplot(111)
plot_orth_proj(abs(k-k1).*maskn,L,L/2,L-100);
title('absolute difference between k and k1')

% 
% figure(5)
% plot_orth_proj(abs(noise).*maskn,L,L/2,L-100);
% title('noise added to g to obtain g1')

figure(7)
plot_orth_proj(mask.*maskn,L,L/2,L-100);
title('Mask applied on g to obtain g1')


