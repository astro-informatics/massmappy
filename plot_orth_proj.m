function [] = plot_orth_proj( f, L, theta0, phi0, plane)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% centered in L/2, L :
% x = -0.57 : 0.001: -0.03;
% y = 0.68 : 0.001: 0.87;

% centered in L/2, L-100 :
x = -0.24 : 0.0005: 0.28;
y = 0.68: 0.0005: 0.87;

T1 = size(x,2);
T2 = size(y,2);

% phi = 0:(2*L-2) *2*pi/(2*L-1);
% theta = 0:(L-1) * pi / L;
% 
% phi1 = phi0 - phi;
% theta1 = theta0 - theta;

fp = zeros(T2,T1);

dec = zeros(1,T2);
ra = zeros(T2,T1);

amin = T1;
amax = 0;
cmin = T2;
cmax = 0;

for k = 1:T2
 
   theta1 = asin(y(k));
   th = round( theta0 + theta1 * L / pi) +1; 
   
   dec(k) = 90 - (th * 180 / L );

   for l = 1:T1
    
   if( abs(x(l)/cos(theta1)) >= 1)
       fp(k,l) = NaN;
   else
       phi1 = asin( x(l)/cos(theta1) );
       ph = round(  phi0 + phi1 *(2*L-1)/(2*pi) );
       ph = mod(ph, 2*L-1)+1;

       ra(k,l) = ph * 180 /(2*L-1);
       
       fp(k,l) = f(th, ph);
       
       if(l < amin)
           amin = l;
       end
       if(l > amax)
           amax = l;
       end
       if(k < cmin)
           cmin = k ;
       end
       if(k > cmax)
           cmax = k;
       end
   end
   end
end
a = amin : amax;
c = cmin : cmax;

if( plane == 1)
    fp1 = nan2zero(fp(c,a));

    kp1 = shear2conv_planorth(fp1,1:size(fp1,1),1:size(fp1,2),15/(3*180)*pi,dec,ra);


    b = imagesc(ra(T2,:), dec, real(-kp1));
    set(gca, 'YDir', 'normal');
    set(b,'AlphaData',~isnan(fp(c,a)));
    
else

    b = imagesc(ra(T2,:), dec, real(fp(c,a)));
    set(gca, 'YDir', 'normal');
    set(b,'AlphaData',~isnan(fp(c,a)));  
end

    xlabel('RA')
    ylabel('Dec')
    
    for k = -45 : -5 : -60
      x = 55:1:90;
      y = k*ones(size(x));
      line(x,y, 'Color', [0.7,0.7,0.7]);
    end
    
    for k = 65 : 5 : 80
      x = 55:0.1:90;
      
      ind = find_closest2(ra(1,:),k);
      
      x1 = ra(T2,ind);
      x2 = k;
      
      y = (x - x2)/(x1-x2) * (dec(1)-dec(T2)) +dec(T2);
      
      line(x,y, 'Color', [0.7,0.7,0.7]);
    end
    
    k = 85;
    x = 55:0.1:90;
    ind = find_closest2(dec,k);
      
    x1 = ra(T2,T1);
    x2 = k;
      
    y = (x - x2)/(x1-x2) * (-49-dec(T2)) +dec(T2);
      
    line(x,y, 'Color', [0.7,0.7,0.7]);    
end

