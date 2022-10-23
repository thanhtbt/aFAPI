function [theta] = DOA_RPAST(X,beta,r)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% r     : target rank (number of signals)

[n, N] = size(X);
theta  = [];
%Initialization
W   = eye(n,r);
P   = eye(r);


N_w    = 14;
beta_s = 0.99;
beta_u = 0.99;
A_s    = zeros(N_w,1);
A_u    = zeros(N_w,1);
sgm_tld= 1;
uu     = 0;
Ti     = 1.96*sqrt(sgm_tld);
cpt    = 1;

%Processing
for k = 1 : N
    y     = W'*X(:,k);
    h     = P*y;
    g     = h/(beta + y'*h);
    e     = X(:,k) - W*y;
    De    = abs(norm(e,'fro') - uu);
    if De < Ti
        P     = P/beta - (g*h')/beta;
        W     = W + e*g';
        Ut{1,k} = W;
    else
        cpt    = cpt + 1;
        Ut{1,k} = W;
        
    end
        
    A_u(1:N_w-1) = A_u(2:N_w);
    A_u(N_w)     = norm(e,'fro');
    uu           = beta_u*uu + (1-beta_u)*median(A_u);
    
    A_s(1:N_w-1) = A_s(2:N_w);
    A_s(N_w)     = (De)^2;
    sgm_tld      = beta_s*sgm_tld + 1.483*(1+ 5/(N_w-1))*(1-beta_s)*median(A_s);
    Ti           = 1.96*sqrt(sgm_tld);
    
    
    % DOA tracking using ESPRIT
    % ESPIRIT algorithm
    D1 = W(1:n-1,:);
    D2 = W(2:n,:);
    
    % TLS solution
    [~, ~, U] = svd([D1 D2]);
    U12 = U(1:r,(r+1):(2*r));
    U22 = U((r+1):(2*r),(r+1):(2*r));
    Phi = -U12*pinv(U22); % TLS solution for Psi
    
    [~,Dest] = eigs(Phi);  
    west = angle((diag(Dest)))/(2*pi);
    
    theta = [theta sort(west)];
    
    
end


