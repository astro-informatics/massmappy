function [ out_ind] = nearest_ind( f, f1 )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T = max(size(f));
out_ind = zeros(size(f));

for k = 1:T
   
    n = 1;
    
    while(f(k) > f1(n) && n < max(size(f1)))
        n = n+1;
    end
    
    if( n > 1 && ( abs(f(k)-f1(n)) < abs(f(k)-f1(n-1)) ) )
        n = n-1;
    end
        
    out_ind(k) = n;
end

end

