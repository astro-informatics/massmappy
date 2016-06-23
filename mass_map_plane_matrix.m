
addpath('/unsafe1/opt/ssht/src/matlab')


rng(2)

N = 100;
L = 2048;
T1 = L * (2*L - 1);
N1 = L/2-N/2+1  : L/2+N/2;
N2 =  L-N/2+1 : L+N/2 ;

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

k = real( ssht_inverse(smooth_lm(klm,pi/(180*3),L),L) );

gamma = ssht_inverse( smooth_lm(gammalm,pi/(3*180),L),L,'Spin',2);
gamma1 = gamma / (max(max(abs(gamma))) - min(min(abs(gamma))));
k = shear2conv(gamma1,L);

M = zeros(size(gamma1));
M(L/2-N/2+1  : L/2+N/2 , L-N/2+1 : L+N/2 ) = ones(N,N);

noise = normrnd(0,1e-02,size(gamma1)) + 1i * normrnd(0,1e-02,size(gamma1)) ;

gamma1 = gamma1 .* M ; %+ noise;
% 
% k1 = shear2conv(gamma,L);





% gamma1 = (gamma +noise) .* mask;

%------------------------- CALCUL MATRICES --------------------------------
% N=1;
%
% C = power_spectrum_computation(Cell, L);
% k2 = mldivide( Om + inv(C) , reshape(Phit( 1/N * gamma),T1,1));
% k2 = reshape(k2, L, 2*L-1);
%
% % k2lm(1:4) = zeros(4,1);
% % 
% % k2 = ssht_inverse(k2lm,L,'Reality',true);
%
% figure(1)
% subplot(211)
% ssht_plot_mollweide(k,L, 'ColourBar', true);
% title('Simulated convergence k')
% 
% subplot(212)
% ssht_plot_mollweide(k2,L, 'ColourBar', true);
% title('Convergence k1 obtained from matrix computation')
% 
% 
% figure(3)
% ssht_plot_mollweide(abs(k-k2),L, 'ColourBar', true);
% title('Difference between k and k1')


% ------------------------- CALCUL PLAN/SPHERE -----------------------------

k2p = shear2conv_plan(gamma1, L/2-N+1 : L/2+N , L-N+1 : L+N, 0, L);
% k21p = shear2conv_plan(gamma,  L/2-N+1 : L/2+N ,  L/2-N+1 : L/2+N  , 20*pi/(180*3), L);

k2p = k2p(N/2+1: 3*N/2, N/2+1: 3*N/2);
k1 = smooth_sph(k, pi/(180*3), L);
k1p = k1(N1 ,N2);
kp = k(N1 , N2);
% k1p = k1(L/2-N/2+1 : L/2+N/2 , L-N/2+1 : L+N/2 );
gp =  real( gamma1(N1 , N2 ) );
% k2p = rand(size(kp)) * 0.04;
% 
sqrt( sum(sum( (real(kp)-real(k2p)).^2))/(N*N))

figure(3)
subplot(221)
imagesc(real(kp))
title('Simulated convergence');
subplot(222)
imagesc(gp)
title('corresponding shear');
subplot(223)
imagesc(real(k2p))
title('plane convergence k2');
subplot(224)
imagesc(abs(real(kp)-real(k2p)))
title('k-k2');

figure(4)
subplot(121)
imagesc(real(kp))
title('Simulated convergence');
subplot(122)
imagesc(real(k2p))
title('plane convergence k2');
% 
% figure(3)
% subplot(121)
% imagesc(abs(kp-k1p))
% title('k-k1');
% subplot(122)
% imagesc(abs(kp-real(k2p)))
% title('k-k2');