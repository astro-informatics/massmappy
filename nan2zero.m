function [ out_f ] = nan2zero( f )

out_f = f;
[t1,t2] = size(f);

for k = 1:t1
    for l = 1:t2
        
        if(isnan(f(k,l)))
            
            out_f(k,l) = 0;
            
        end
    end
end


end

