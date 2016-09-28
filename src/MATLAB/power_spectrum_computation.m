function [ out_C ] = power_spectrum_computation( Cell, L )

% Computes power spectrum of kappa.

T1 = (2*L-1) * L;
x = zeros(T1,T1);
out_C = zeros(T1,T1);

m = 0;
n=0;
    
    for thm = (0: 1/(L-1) : 1)*pi  
    for phm = (0: 1/(2*L-1) : 1-1/(2*L-1))*2*pi    
        
        m = m+1
        for thn = (0: 1/(L-1) : 1)*pi  
        for phn = (0: 1/(2*L-1) : 1-1/(2*L-1))*2*pi 
            
            n=n+1;
            x(m,n) = cos(thn)*cos(thm) + sin(thn)*sin(thm)*(cos(phm)*cos(phn)+sin(phm)*cos(phn));
        end
        end
        n=0;
    end
    end
    
for l = 1:L
    l
    out_C = out_C + Cell(l) * legendreP(l,x);
    
end

end

