function [ out_noise ] = noise_per_pixel( k, L )

%   k is convergence, and L the number of samples.

Ngal = 1e09;
mu = 1;

dtheta = pi /(L-1);
dphi = 2*pi /(2*L-1);
dOmega = dtheta * dphi / (4* pi);

theta = dtheta * (0 : 1 : L-1);

Dtheta = diag(sin(theta));

% Nk est la norme 1 de k surla sphere. On doit donc multiplier par
% sin(theta) toutes les lignes de k (qui dépendent de theta).

Nk = sum(sum(Dtheta * (mu + k)))*dtheta*dphi;

lambda = (mu+k) * dOmega * Ngal / Nk;

n = poissrnd(lambda);

sum(sum(n))

%out_noise = 0.2 * sqrt(1 ./ n);
sigma = 0.2 * sqrt(1 ./ n);

out_noise = 1/sqrt(2) * (normrnd(zeros(size(sigma)),sigma) + i*normrnd(zeros(size(sigma)),sigma));
end

