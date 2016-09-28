function [alpha, beta, gamma ] = angle2notrhpole( l, k, L )

    theta = ( 0:(L-1)) *pi/L;
    phi = ( 0:(2*L-2) ) * 2*pi/ (2*L-1);
 
    gamma = -phi(k);
    alpha = pi/2;
    beta = pi-theta(l);

end

