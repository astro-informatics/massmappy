addpath('/unsafe1/opt/ssht/src/matlab')
load('DES1')
L = 512;
g = zeros(L,2*L-1);
w = zeros(L,2*L-1);
nbGal = zeros(L,2*L-1);

ra1 = RA / 180 * (2*L-1);
dec1 = (90 - dec) /180 * L;

T = size(e1,1);

amin = L;
amax = 0;
cmin = 2*L-1;
cmax = 0;

for l = 1:T
    g( floor(dec1(l)), floor(ra1(l))) = g( floor(dec1(l)), floor(ra1(l))) + weight(l)*( e1(l)-c1(l) + 1i * (e2(l)-c2(l) ) );
    w( floor(dec1(l)), floor(ra1(l))) = w( floor(dec1(l)), floor(ra1(l))) + weight(l)*(1+mcorr(l));
    nbGal(floor(dec1(l)), floor(ra1(l))) = nbGal( floor(dec1(l)), floor(ra1(l)))  + 1;
    
    if(floor(dec1(l)) < amin)
        amin = floor(dec1(l));
    end
    if(floor(dec1(l)) > amax)
        amax = floor(dec1(l));
    end
    if(floor(ra1(l)) < cmin)
        cmin = floor(ra1(l));
    end
    if(floor(ra1(l)) > cmax)
        cmax = floor(ra1(l));
    end
end
    
M = double(w == 0);
w = w + M;

m1 = 0;
m2 = 0;

for l1 = 1:L
    for l2 = 1:2*L-1
        if(M(l1,l2) == 0)
            m1 = m1 + l1;
            m2 = m2 + l2;
        end
    end
end

m1 = round(m1/sum(sum(1-M)));
m2 = round(m2/sum(sum(1-M)));

g = g ./ w;


% g = g .* double( abs(g) < 0.2);
% 
% mask = double(nbGal > 100);
mask = 1-double(nbGal == 0);

g = g .* mask;

figure(3)
ssht_plot_mollweide(g.*mask2nan(mask),L,'ColourBar', false);
% a = 189:215;
% c = 170:260;

maskn = mask2nan(mask);
a = amin:amax;
c = cmin:cmax;

% coord matrix :  (185:220,140:300)