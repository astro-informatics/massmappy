function [ out_n ] = rayleigh_distrib(N)

out_n = 0;

A = 1/0.8309;
emax = 0.804;
e0 = 0.0256;
a = 0.2539;

p =@(e) A * e .* (1 - exp( (e-emax)./a))./((1+e).*sqrt((e.^2+e0^2)));

pmax = 1;
    
k = 0;
% 
% while( k < N )
%     
%     a1 = pmax*rand();
%     b1 = emax*rand();
%         
%     if(p(b1) > a1)
%        k = k+1;
%        phi = 2 * pi * 1i * rand();
% 
%        out_n = out_n + b1 .* exp(phi);
%     end
% end
% 
% out_n = out_n/N;

N1 = N;

while(N1 > 0)
    a1 = pmax * rand(1,N1);
    b1 = emax * rand(1,N1);
    M = double((p(b1) > a1) );
    
    phi = exp( 2 * pi * 1i * rand(1,N1));
    
    out_n = out_n + sum( sum( b1 .*M .* phi ))/N;
    
    N1 = N1 - sum(sum(M));
end


end

