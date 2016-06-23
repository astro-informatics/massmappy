function [ out_ind ] = find_closest2( tab, n )

%tab needs to be ordered !
k = 1;

while(tab(k) < n && k < max(size(tab)) )
    
    k = k+1;
    
end
    out_ind = k;

end

