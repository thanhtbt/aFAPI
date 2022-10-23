function [Ut,rho,eta] = alpha_FAPI(X,beta,alpha,U_tr,lp)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% alpha : divergence
% U_tr  : the set of true subspaces with time
% lp    : (0 < lp <=2)

% Author: Le Trung Thanh  (University of Orleans, France)
% Email : thanhle88.tbt@gmail.com
% Cite  : L.T. Thanh, A.M. Rekavandi, S. Abd-Krim, & K. Abed-Meraim.
          ... “Robust Subspace Tracking With Contamination via Alpha-Divergence”. 
          ... ICASSP 2023 (submitted).
    
if nargin <= 4
    lp = 1.5;
else
end

% Initialization
[n, N] = size(X);
r      = min(size(U_tr{1,1}));
W      = eye(n,r);
Z      = eye(r);
c      = (1-alpha)/2;

% Performance Metrics
rho = zeros(1,N);
eta = zeros(1,N);

%% Tracking ...
for k  = 1:N
    %% Main Program
    y   = W'*X(:,k);
    h   = Z*y;
    e   = X(:,k) - W*y;
    w   = exp(-c*norm(e)^lp);
    g   = h*w/(beta + y'*h*w);
    gn  = norm(g)^2;
    eps =  norm(r)^2;  % (norm(X(:,k),'fro'))^2 - (norm(y,'fro'))^2; %
    tau = eps/(1+ eps*gn + sqrt(1 + eps*gn) );
    mu  = 1 - tau*gn;
    yp  = mu*y + tau*g;
    hp  = Z'* yp;
    r   = (tau/mu)*(Z*g - (hp'*g)*g);
    Z   = (1/beta)* (Z - g*hp' + r*g');
    ep  = mu* X(:,k) - W*yp;
    W   = W + ep*g';
    
    %% Performance Evaluation
    V       = orth(U_tr{1,k}); % true subspace at time t
    rho(k)  = abs(trace(W'*(eye(n)-V*V')*W)/trace(W'*(V*V')*W));
    eta(k)  = sin(subspace(W,V));
    Ut{1,k} = W;
end
end

