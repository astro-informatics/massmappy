function [] = plot_gnomonic_proj(f, L, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(nargin > 2)
    alpha = varargin{1};
    beta = varargin{2};
    gamma = varargin{3};
    f = rotation_sphere(f,alpha,beta,gamma,L);
    
    if(nargin > 5)
      mask = varargin{4};       
      mask = rotation_mask(mask,alpha,beta,gamma,L);
      mask = mask2nan(mask);
      f = f .*mask;
    end
end

    
x = -0.3:0.001:0.3;
y = -0.3:0.001:0.3;

T = size(x,2);
% theta_f = 0: 1 : (L-1) * pi/L;
% phi_f = 0:1:(2*L-2) * 2*pi/(2*L-1);

z = ones(T,1)*x + 1i * y' * ones(1,T);
theta = atan(abs(z));
theta = round((pi-theta) * L/pi)+1;

phi = angle(z) + pi;
phi = round( phi * (2*L-1)/(2*pi) )+1;

phi = phi .* (phi <= 2*L-1) + ones(T,T).* (1-(phi <= 2*L-1));
theta = theta .* (theta <= L) + ones(T,T).* (1-(theta <= L));

fp = zeros(T,T);


for k = 1:T
   for l = 1:T
        
       if( abs(z(k,l)) > 1) 
           fp(k,l) = NaN;
       else
        fp(k,l) = f(theta(k,l),phi(k,l));
       end
       
   end
end



    b = imagesc(real(fp));
    set(b,'AlphaData',~isnan(fp));
end

