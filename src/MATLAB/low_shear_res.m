function [ g_l ] = low_shear_res( g_h, L)

%  g_h is high resolution shear. From now we suppose g_h has 4*L samples
%  g_l is low resolution shear with resolution L.

g_l = zeros(L, 2*L-1);
l = floor((1:2*L)*(8*L-1)/(2*L-1));

phi_h = 2*pi/(8*L-1) *(0:8*L+5);
theta_h = pi/(4*L)*(0:4*L);

phi_l = 2*pi/(2*L-1) *(0:2*L-1);
theta_l = pi/L*(0:L);

for j0 = 1:L
    dOmega = pi/L * (sin(theta_l(j0)) + sin(theta_l(j0+1)))/2;
    
    for k0 = 1: 2*L-1

        n = l(k0+1) - l(k0) + 1;
        
        w = ones(1,n);
        
        
        
        w(1) = ( phi_h(l(k0)) - phi_l(k0) )/( phi_h(l(k0)) - phi_h(l(k0)-1) );
        
        w(n) =  ( phi_l(k0+1) - phi_h(l(k0+1)-1) )/( phi_h(l(k0+1)) - phi_h(l(k0+1)-1) );
        
        s = 0;
        
        for j = 4*(j0-1)+1 : 4*(j0-1) + 4
            
            if(l(k0+1)-1 <= 8*L-1);
                s = s + ( sin(theta_h(j+1)) + sin(theta_h(j)) )/2 * sum( w .* g_h( j, l(k0)-1:l(k0+1)-1 ) );
            else
                ind = mod( l(k0+1)-1, 8*L-1);
                m = 8*L-1-(l(k0)-1)+1;
                s = s + ( sin(theta_h(j+1)) + sin(theta_h(j)) )/2 * sum( w(1:m) .* g_h( j, l(k0)-1:8*L-1 ) );
                s = s + ( sin(theta_h(j+1)) + sin(theta_h(j)) )/2 * sum( w(m+1:n) .* g_h( j, 1:ind ) );
            end

        end
        
        g_l(j0,k0) = 2*pi / (2*L-1) * s / dOmega; 
        
    end
    
end

g_l = max(abs(g_h(:)))/max(max(abs(g_l))) * g_l;


end

