
addpath('/unsafe1/opt/s2let/src/main/matlab')
addpath('/unsafe1/opt/ssht/src/matlab')

L = 1024;
B = 3;
J_min = 2;
nside_recon = 128;
s2let_path = '/home/ap2';
inputfile = strcat(s2let_path,'/Euclid_footprint_and_StarMask4096.fits');

% Read the file
[f_ini, nside] = s2let_hpx_read_real_map(inputfile);

% Band limit the data
flm = s2let_hpx_map2alm(f_ini, 'L', L);
f = s2let_hpx_alm2map(flm, nside_recon, 'L', L);

% g1 = ssht_inverse(flm,L, 'Reality', true);

% 
% % Perform decomposition
% [f_wav, f_scal] = s2let_transform_axisym_analysis_hpx(f,'B',B,'L',L,'J_min',J_min);
% 
% % Plot
% J = s2let_jmax(L, B);
% zoomfactor = 1.2;
% ns = ceil(sqrt(2+J-J_min+1)) ;
% ny = ns - 1 ;
% nx = ns ;
% figure('Position',[100 100 1300 1000])
% 
% subplot(nx, ny, 1);
% s2let_hpx_plot_mollweide(f);
% campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
% title('Initial band-limited data')
% 
% subplot(nx, ny, 2);
% s2let_hpx_plot_mollweide(f_scal);
% campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
% title('Scaling fct')
% 
% for j = J_min:J
%    subplot(nx, ny, j-J_min+3);
%    s2let_hpx_plot_mollweide(f_wav{j-J_min+1});
%    campos([0 0 -1]); camup([0 1 0]); zoom(zoomfactor)
%    title(['Wavelet scale : ',int2str(j)-J_min+1])
% end 
