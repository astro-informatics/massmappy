function [error_s, error_p,rms_k] = mass_map_plane_function( Cell, L, N, buffer_ratio,theta0, rect_ratio, band_limit_mask, mask_type)


addpath('/unsafe1/opt/ssht/src/matlab')


rng(4)

N_p = buffer_ratio*N;

N1 = L/2+theta0-N/2+1  : L/2+theta0+N/2;
N2 = L- floor(N*rect_ratio)/2+1 : L+floor(N*rect_ratio)/2;

N1p = L/2+theta0-N_p/2+1  : L/2+theta0+N_p/2;
N2p =   L-floor(N_p*rect_ratio)/2+1 : L+floor(N_p*rect_ratio)/2;

sigma_smooth = 0;%pi/(2*L);

D = zeros(L,1);

klm = [];
gammalm= [];

for l = 1:L
    D(l+1) = -1 /(l*(l+1)) * sqrt((l+2)*(l+1)*l*(l-1));
end

for l = 0:(L-1)  
    
    km = zeros(2*l+1,1);
    km = normrnd(0, sqrt( Cell(l+1) ), [ l ,1]) + 1i * normrnd(0, sqrt( Cell(l+1) ), [ l ,1]) ;
    km(l+1,1) = normrnd(0, sqrt( Cell(l+1) ));

    for m = 1:l
        km(l+1+m,1) = (-1).^(m) * conj(km(l+1-m,1)); 
    end
%    km = normrnd(0, sqrt( Cell(l+1) ), [ 2*l+1 ,1]) + 1i * normrnd(0, sqrt( Cell(l+1) ), [ 2*l+1 ,1]) ;

    klm = [klm; km ];
    
    gammalm = [gammalm; D(l+1)*km];
    
end

klm(1:4) = [0;0;0;0];

%k = real( ssht_inverse(smooth_lm(klm,sigma_smooth,L),L) );
%figure, ssht_plot_mollweide(k, L, 'Method', 'MW','Mode',0);
%drawnow;

%gamma = ssht_inverse( smooth_lm(gammalm,pi/(3*180),L),L,'Spin',2);

k = real( ssht_inverse(klm,L) );



if band_limit_mask == 1,
    M = zeros(size(k));
    if mask_type == 0,    
        M(N1 , N2 ) = ones(N,floor(N*rect_ratio));
    else
        ang = pi*N/(2*L);
        [thetas, phis] = ssht_sampling(L, 'Grid', true);
        d = acos(sin(thetas).*sin(thetas(L/2+1,L+1)).*cos(phis-phis(L/2+1,L+1))+cos(thetas).*cos(thetas(L/2+1,L+1)));
        M(d<ang) = 1.0;
    end
    k = k.*M;
    k_old = k;
    k_lm = ssht_forward(k,L,'Reality',true);
    for el=0:(L-1) 
        
        gammalm = [gammalm; D(l+1)*km];
    %M_lm = recsph_gaussian_smooth(M_lm, L, 2*pi/L);
    k = ssht_inverse(k_lm,L,'Reality',true);
    figure(10)
    ssht_plot_mollweide(k,L)
    drawnow;
    figure(11)
    ssht_plot_mollweide(k-k_old,L)
    drawnow;
end


gamma = ssht_inverse(gammalm,L,'Spin',2);

gamma1 = gamma;
% gamma1 = (real(gamma) - min(min(real(gamma))))/ (max(max(real(gamma))) - min(min(real(gamma)))) - 0.5;
% gamma1 = gamma1 + 1i * ((imag(gamma)- min(min(imag(gamma))))/ (max(max(imag(gamma))) - min(min(imag(gamma)))) - 0.5);
% k = shear2conv(gamma1,L);
%min(min(real(gamma)))
% max(max(real(gamma)))
% 

if (band_limit_mask == 0),
    M = zeros(size(gamma1));
    if mask_type == 0,    
        M(N1 , N2 ) = ones(N,floor(N*rect_ratio));
    else
        ang = pi*N/(2*L);
        [thetas, phis] = ssht_sampling(L, 'Grid', true);
        d = acos(sin(thetas).*sin(thetas(L/2+1,L+1)).*cos(phis-phis(L/2+1,L+1))+cos(thetas).*cos(thetas(L/2+1,L+1)));
        M(d<ang) = 1.0;
    end
    gamma1 = gamma1 .* M ;
end
% 
if band_limit_mask == 2,
    M = zeros(size(gamma1));
    if mask_type == 0,    
        M(N1 , N2 ) = ones(N,floor(N*rect_ratio));
    else
        ang = pi*N/(2*L);
        [thetas, phis] = ssht_sampling(L, 'Grid', true);
        d = acos(sin(thetas).*sin(thetas(L/2+1,L+1)).*cos(phis-phis(L/2+1,L+1))+cos(thetas).*cos(thetas(L/2+1,L+1)));
        M(d<ang) = 1.0;
    end
    M_old = M;
    M_lm = ssht_forward(M,L,'Reality',true);
    M_lm = recsph_gaussian_smooth(M_lm, L, 50*pi/L);
    M = ssht_inverse(M_lm,L,'Reality',true);
    figure(10)
    ssht_plot_mollweide(M,L)
    drawnow;
    figure(11)
    ssht_plot_mollweide(M-M_old,L)
    drawnow;
    gamma1 = gamma1 .* M ;
end


%noise = normrnd(0,1e-02,size(gamma1)) + 1i * normrnd(0,1e-02,size(gamma1)) ;
%gamma1 = gamma1 + noise;
% figure, ssht_plot_mollweide(real(gamma1), L, 'Method', 'MW','Mode',0);

k1 = shear2conv_new(gamma1,L);

% ------------------------- CALCUL PLAN/SPHERE -----------------------------
%k2p = shear2conv_plan(gamma1, L/2-N+1 : L/2+N , L-N+1 : L+N, 10*sigma_smooth, L); 
[ k2p] = shear2conv_plan(gamma1, N1p , N2p, 0.0, L); %10*sigma_smooth
% k21p = shear2conv_plan(gamma,  L/2-N+1 : L/2+N ,  L/2-N+1 : L/2+N  , 20*sigma_smooth, L);


k2p = k2p(N_p/2-N/2+1: N_p/2+N/2, (N_p/2-N/2)*rect_ratio+1: (N_p/2+N/2)*rect_ratio);

% k2p_s = zeros(L,2*L-1);
% k2p_s(N1,N2) = k2p;
% k2p_s = smooth_sph(k2p_s, sigma_smooth, L);
% k2p =  k2p_s(N1,N2);

ks = smooth_sph(k, sigma_smooth, L);
%ks = smooth_sph(k, 0.0, L);
ksp = ks(N1 ,N2);

kp = k(N1 , N2).*M(N1 , N2);

k1s = smooth_sph(k1, sigma_smooth, L);
%k1s = smooth_sph(k1, 0.0, L);%sigma_smooth
k1sp = k1s(N1 ,N2);
k1p = k1(N1,N2);

% k1p = k1(L/2-N/2+1 : L/2+N/2 , L-N/2+1 : L+N/2 );
gp = smooth_sph_gamma(gamma, sigma_smooth, L);
gp =  real( gp(N1 , N2 ) );

% k2p = rand(size(kp)) * 0.04;
% 
rms(real(kp)-real(k2p));

figure(3)
subplot(321)
imagesc(real(kp))
title('Simulated convergence');
subplot(322)
imagesc(gp)
title('corresponding shear');
subplot(323)
imagesc(real(k1p))
title('sherical convergence k1');
subplot(324)
imagesc(abs(real(kp)-real(k1p)))
title('k-k1');
subplot(325)
imagesc(real(k2p))
title('plane convergence k2');
subplot(326)
imagesc(abs(real(kp)-real(k2p)))
title('k-k2');


figure(4)
subplot(121)
imagesc(real(ksp))
title('Simulated convergence');
subplot(122)
imagesc(real(k2p))
title('plane convergence k2');

drawnow;
error_s =  rms(real(kp(:))-real(k1p(:)));
error_p =  rms(real(kp(:))-real(k2p(:)));
rms_k = rms(real(kp(:)));
% 
% figure(3)
% subplot(121)
% imagesc(abs(kp-k1p))
% title('k-k1');
% subplot(122)
% imagesc(abs(kp-real(k2p)))
% title('k-k2');
