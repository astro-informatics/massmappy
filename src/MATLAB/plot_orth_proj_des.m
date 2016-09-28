function [] = plot_orth_proj_des()

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load('DES.mat')

% centered in L/2, L :
% x = -0.57 : 0.001: -0.03;
% y = 0.68 : 0.001: 0.87;

% centered in L/2, L-100 :
% x = -0.28 : 0.0025: 0.22;
% y = 0.47: 0.0025: 0.73;
% x = -1:0.01:1
% y = -1:0.01:1
x = -0.19 : 0.001: 0.07;
y = -0.88: 0.001: -0.68;
 
T1 = size(x,2);
T2 = size(y,2);

% phi = 0:(2*L-2) *2*pi/(2*L-1);
% theta = 0:(L-1) * pi / L;
% 
% phi1 = phi0 - phi;
% theta1 = theta0 - theta;

fp = zeros(T2,T1);

dec1 = zeros(1,T2);
ra1 = zeros(T2,T1);

amin = T1;
amax = 0;
cmin = T2;
cmax = 0;

for k = 1:T2
 
   theta1 = acos(y(k));
%    th = round( theta0 + theta1 * L / pi) +1; 
   
   dec1(k) = theta1 * 180 / pi - 90;

   for l = 1:T1
    
   if( sin(theta1) == 0 || abs(x(l)/sin(theta1)) >= 1)
       fp(k,l) = NaN;
   else
       phi1 = asin( x(l)/sin(theta1) );
%        ph = round(  phi0 + phi1 *(2*L-1)/(2*pi) );
%        ph = mod(ph, 2*L-1)+1;

       ra1(k,l) = phi1 * 180 /pi + 80;

   end
   end
end

T = max(size(e1));

% 
% k1 = nearest_ind(abs(dec),abs(dec1));
% 
% k2 = zeros(1,T);
% for l = 1:T
%     RA(l);
%     k2(l) = nearest_ind(RA(l), ra1(k1(l),:));
%     ra1(k1(l),k2(l));
% end

g = zeros(T2,T1);
w = zeros(T2,T1);
nbGal = zeros(T2,T1);

for l = 1:T
    dec(l);
    th = (90-dec(l))*pi/180;
    ph = (RA(l)+100)*pi/180;
    y1 = cos(th);
    x1 = -sin(th)*sin(ph);
    k2 = nearest_ind(x1,x);
    k1 = nearest_ind(y1,y);
    dec1(k1);
    g( k1, k2) = g( k1, k2) + weight(l)*( e1(l)-c1(l) + 1i * (e2(l)-c2(l)) );
    w( k1, k2) = w( k1, k2) + weight(l)*(1+mcorr(l));
    nbGal( k1, k2) = nbGal( k1, k2)  + 1;
    
end

M = double(w == 0);
w = w + M;
g = g ./ w;

max(max(abs(g)))
mask = (1-double(nbGal == 0)).* (1-double( isnan(fp)));
% mask = (1-double( isnan(fp)));
max(max(-mask))
max(k1);
max(k2);
maskn = mask2nan(mask);

a = amin : amax;
c = cmin : cmax;

abs(ra1(1,1)-ra1(1,T1))
abs(ra1(T2,1)-ra1(T2,T1))

    kp1 = shear2conv_planorth(g,1:T2,1:T1,15/(3*180)*pi,dec1,ra1);


    b = imagesc(ra1(1,:), -dec1, -real(kp1).*maskn);

%     b = imagesc(x, y, real(kp1));    
    
 %---------------------DETAILS D'AFFICHAGE----------------------------
    
    
    set(gca, 'YDir', 'normal');
    set(b,'AlphaData',~isnan(maskn));
 

    
    xlabel('RA')
    ylabel('Dec')
    
    
%     
%     for k = -45 : -5 : -60
%       x = 55:1:90;
%       y = k*ones(size(x));
%       line(x,y, 'Color', [0.7,0.7,0.7]);
%     end
%     
%     for k = 65 : 5 : 80
%       x = 55:0.1:90;
%       
%       ind = find_closest2(ra1(1,:),k);
%       
%       x1 = ra1(T2,ind);
%       x2 = k;
%       
%       y = (x - x2)/(x1-x2) * (dec1(1)-dec1(T2)) +dec1(T2);
%       
%       line(x,y, 'Color', [0.7,0.7,0.7]);
%     end
%     
%     k = 85;
%     x = 55:0.1:90;
%     ind = find_closest2(dec1,k);
%       
%     x1 = ra1(T2,T1);
%     x2 = k;
%       
%     y = (x - x2)/(x1-x2) * (-49-dec1(T2)) +dec1(T2);
%       
%     line(x,y, 'Color', [0.7,0.7,0.7]);    
end

