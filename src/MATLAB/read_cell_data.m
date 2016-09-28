addpath('/unsafe1/opt/ssht/src/matlab')

fichier = fopen('cell_icosmo.txt','r');

L = 2048;

C = textscan(fichier,'%f64 %f64 %f64','EmptyValue',0);

l1=C{1};
C1 = C{2};
Cell = zeros(1,L);

j = 1;
t = 1;

for l = 1:L
    
    while( l1(j) < l)
        
        j = j+1;
    end
    
    
    if(l1(j) == l)
        Cell(j) = C1(j);
    else
        t = (l-l1(j) )/(l1(j-1)-l1(j));
        Cell(l) = t*C1(j-1) + (1-t)*C1(j);
    end
    
end

