function [Ut,rho,eta] = RPAST_Tracking(X,beta,U_tr)
% X    : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta : forgetting factor
% U_tr : the set of true subspaces with time


[n, N] = size(X);
r      = min(size(U_tr{1,1}));

%Initialization
W   = eye(n,r);
P   = eye(r);
rho = zeros(1,N);
eta = zeros(1,N);

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
    V     = orth(U_tr{1,k}); % true subspace at time t
    y     = W'*X(:,k);
    h     = P*y;
    g     = h/(beta + y'*h);
    e     = X(:,k) - W*y;
    De    = abs(norm(e,'fro') - uu);
    if De < Ti
        P     = P/beta - (g*h')/beta;
        W     = W + e*g';
        rho(k)  = abs(trace(W'*(eye(n)-V*V')*W)/trace(W'*(V*V')*W));
        eta(k)  = sin(subspace(W,V));
        Ut{1,k} = W;
    else
        cpt    = cpt + 1;
        rho(k)  = abs(trace(W'*(eye(n)-V*V')*W)/trace(W'*(V*V')*W));
        eta(k)  = sin(subspace(W,V));
        Ut{1,k} = W;
        
    end
    
    
    A_u(1:N_w-1) = A_u(2:N_w);
    A_u(N_w)     = norm(e,'fro');
    uu           = beta_u*uu + (1-beta_u)*median(A_u);
    
    A_s(1:N_w-1) = A_s(2:N_w);
    A_s(N_w)     = (De)^2;
    sgm_tld      = beta_s*sgm_tld + 1.483*(1+ 5/(N_w-1))*(1-beta_s)*median(A_s);
    Ti           = 1.96*sqrt(sgm_tld);
    
end


