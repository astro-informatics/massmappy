
figure(2)
addpath('/unsafe1/opt/ssht/src/matlab')

fichier = fopen('EuclidMask9.txt','r');

L = 2048;

C = textscan(fichier,'%f64 %f64','EmptyValue',0);
T=C{1} + C{2};

T1 = zeros(1,ssht_elm2ind(L-1,L-1));

n = 0;

for l = 0 : L-1
    
    n = n+1;
    
    T1(ssht_elm2ind(l,0)) = T(n);
    
    for m = 1 : l
        
        n = n+1;
        T1(ssht_elm2ind(l,m)) = T(n);
        T1(ssht_elm2ind(l,-m)) = (-1)^m * conj(T(n)) ;
    end
    
end

fclose(fichier);

M = ssht_inverse(T1,L,'Reality',true);
ssht_plot_mollweide(M,L);



