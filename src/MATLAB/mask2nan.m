function [ m ] = mask2nan( m )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[t1,t2] = size(m);

for l = 1:t1
    for k = 1:t2
        
        if( m(l,k) == 0)
            m(l,k) = NaN;
        end
    end
end

end

