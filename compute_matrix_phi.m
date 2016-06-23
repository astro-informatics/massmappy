
addpath('/unsafe1/opt/ssht/src/matlab')

L = 32;
T1 = L * (2*L-1);
T2 = ssht_elm2ind(L-1,L-1);


Gamma = @(x) ssht_forward(x,L); %prend une matrice et calcule sa ssht

Gammat = @(x) ssht_forward_adjoint(x,L); 

M = eye(T1) ; %quick_mask(mask);

d = vector_D(L);

Lambda0 = @(x) ssht_inverse(x,L);

Lambda0t = @(x) ssht_inverse_adjoint(x,L);

Lambda = @(x) ssht_inverse(x,L,'Spin',2);

Lambdat = @(x) ssht_inverse_adjoint(x,L, 'Spin', 2);

N = 1;

% Phi= @(x) M * reshape( Lambda ( d .* Gamma(x) ), T1,1);
% Phit = @(x) Gammat(d .* Lambdat (reshape(M'* x,L,2*L-1)));

Phi= @(x)  Lambda( d.* Gamma (x)) ;
Phit = @(x) Gammat( d.*Lambdat(x)) ;

Om = zeros(T1,T1);

v1 = zeros(L,2*L-1);
ind = 1;
v1(1,1) = 1;
v = Phit(N*Phi(v1));
Om(:,ind) = v(:);

for ind = 2:T1
    
    v1( floor((ind-2)/(2*L-1)) + 1, mod(ind-2,2*L-1)+1) = 0;
    
    v1( floor((ind-1)/(2*L-1)) + 1, mod(ind-1,2*L-1)+1) = 1;
    
    
    v = Phit(N*Phi(v1));
    Om(:,ind) = v(:);
    
    if(mod(ind, 50) == 0)
       ind
    end
end

% 
% Om = zeros(T2,T2);
% v1 = zeros(T2,1);
% 
% ind = 1;
% v1(1) = 1;
% v =  Phit(N*Phi(v1));
% Om(:,ind) = v(:);
% 
% for ind = 2:T2
%     
%     v1(ind-1) = 0;
%     
%     v1(ind) = 1;
%     
%     
%     v = Phit(N*Phi(v1));
%     Om(:,ind) = v(:);
%     
%     if(mod(ind, 50) == 0)
%        ind
%     end
% end
