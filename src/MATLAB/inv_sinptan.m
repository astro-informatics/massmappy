function [ out_inv ] = inv_sinptan(x)

 l = pi * (-0.5 : 0.001 : 0.5);
 T = size(l,2);
 fl = 1./cos(l) + tan(l);
 
 a = 1;
 b = T;
 
 while( b-a > 1)
     
     if(x > fl(floor((a+b)/2)) )
        a = floor((a+b)/2);
     else
        b = floor((a+b)/2);
     end
 end
    
    t = (fl(b)-x)/(fl(b)-fl(a));
    out_inv = t * l(a) + (1-t) * l(b);

end

