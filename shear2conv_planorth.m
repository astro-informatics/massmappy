function [ out_kp ] = shear2conv_planorth( g, a,b, sigma, dec, ra )

% shear2conv_plan converts shear gamma considered on a NxN plane to
% convergence.

% if( mod(N,2) == 1 ) % On va centrer le plot au milieu de la sphere donc N doit être pair.
%     N = N-1;
% end



M = size(a,2);
N = size(b,2);

angM = abs(dec(M)-dec(1));

gp = g(a , b );

%-------------------------------------------------------------------------------------------------------

D = zeros(M, N);

for l1 = 1 : M
        
%   theta1 = (l1 + a(1)-1 )*pi/L;

    angN = abs(ra(l1,N)-ra(l1,1)); 
    
    for l2 = 1 : N
        
%       theta2 = (l2 )*2*pi/(2*L-1);

        l1n = (l1-1-M/2)/angM*M ;
        l2n = (l2-1-N/2)/angN*N ;
        
        
%         D(l1,l2) = (l1n^2 - l2n^2 + 2* 1i *l1n*l2n)/(l1n^2+l2n^2);
        if ~(abs(l1n) < 1E-3 && abs(l2n) < 1E-3)
            D(l1,l2) = (l1n^2 - l2n^2 + 2* 1i *l1n*l2n)/(l1n^2+l2n^2);
        end
        
%       D(l1,l2) = (l1n^2+l2n^2)/(2*l1n*l2n);

    end
end

%-------------------------------------------------------------------------------------------------------

l1 = ones(N,1) * (-M/2+1:M/2) ;
l2 = (-N/2+1:N/2)' * ones(1,M);

smoo = exp( - sigma^2 * ( l1'.^2 + l2'.^2) );

fgp = fftshift( fft2(gp) );

fkp =  conj(D)  .* fgp .* smoo;

out_kp = ifft2(fftshift(fkp) );


end

