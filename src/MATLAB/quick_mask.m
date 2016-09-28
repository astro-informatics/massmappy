function [ out_M ] = quick_mask( m )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

m = m(:);

T = size(m,1);

out_M = zeros(T,T);

n = 0;
v =  zeros(1,T);

for ind = 1:T 
    
    if(ind > 1)
        v(ind-1) = 0;
    end
    
    if m(ind) == 1
       n = n+1;
       v(ind) = 1;
       out_M(n,:) = v;
    end
end

out_M = out_M(1:n,:);

end

