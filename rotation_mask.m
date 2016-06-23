function [ out_f ] = rotation_mask( f, alpha, beta, gamma, L )

% Define parameters.
type = 'colour';

% Generate spherical harmonics.
flm = ssht_forward(f,L, 'Reality', true);

% Compute sampling grids.
%[thetas, phis, n, ntheta, nphi] = ssht_sampling(L, 'Grid', true);

% Precompute Wigner small-d functions
d = zeros(L, 2*L-1, 2*L-1);
d(1,:,:) = ssht_dl(squeeze(d(1,:,:)), L, 0, beta);
for el = 1:L-1
    d(el+1,:,:) = ssht_dl(squeeze(d(el,:,:)), L, el, beta);
end

% Rotate spherical harmonic
flm_rot = ssht_rotate_flm(flm, d, alpha, gamma);

out_f = ssht_inverse(flm_rot,L, 'Reality', true);

out_f= double(out_f > 0.7);


end

