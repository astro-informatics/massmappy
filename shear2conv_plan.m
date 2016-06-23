function [ out_kp ] = shear2conv_plan( g, a,b, sigma, L )

% shear2conv_plan converts shear gamma considered on a NxN plane to
% convergence.

% if( mod(N,2) == 1 ) % On va centrer le plot au milieu de la sphere donc N doit être pair.
%     N = N-1;
% end

M = size(a,2);
N = size(b,2);
Max = max(M,N);

gp = g(a , b );

%-------------------------------------------------------------------------------------------------------

D = zeros(M, N);

for l1 = 1 : M
        
    theta1 = (l1 + a(1)-1 )*pi/L;
    
    for l2 = 1 : N
        
        theta2 = (l2 )*2*pi/(2*L-1);
        
        l1n = l1;
        l2n = l2;
        
%       D(l1,l2) = (l1^2 - l2^2 + 2* 1i *l1*l2)/(l1^2+l2^2);
        
        D(l1,l2) = (l1n^2 - l2n^2 + 2* 1i *l1n*l2n)/(l1n^2+l2n^2);
%       D(l1,l2) = (theta1-1i*theta2)^2/(theta1^2+theta2^2);
    end
end

%-------------------------------------------------------------------------------------------------------

l1 = ones(N,1) * (-M/2+1:M/2) ;
l2 = (-N/2+1:N/2)' * ones(1,M);

smoo = exp( - sigma^2 * ( l1'.^2 + l2'.^2) );

fgp = fft2(gp);

fkp = conj(D) .* fgp .* fftshift(smoo);

out_kp = ifft2( fkp);


end

