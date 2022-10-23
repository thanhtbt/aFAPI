function [theta] = DOA_alpha_FAPI(X,beta,r,alpha,lp)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% r     : target rank (number of signals)
% alpha : divergence
% lp    : (0 < lp <=2 )

% Author: Le Trung Thanh  (University of Orleans, France)
% Email : thanhle88.tbt@gmail.com
% Cite  : L.T. Thanh, A.M. Rekavandi, S. Abd-Krim, & K. Abed-Meraim.
          ... “Robust Subspace Tracking With Contamination via Alpha-Divergence”. 
          ... ICASSP 2023 (submitted).
    
if nargin <= 4
    lp     = 1.5;
else
end

[n, N] = size(X);
W      = eye(n,r);
Z      = eye(r);
P      = eye(n);
c      = (1-alpha)/2;

theta = [];

for k = 1:N
    
    %% aFAPI Section
    y   = W'*X(:,k);
    h   = Z*y;
    w   = exp(-c*norm(X(:,k)-W*y)^lp);
    g   = h*w/(beta + y'*h*w);
    gn  = norm(g)^2;
    eps = norm(X(:,k)-W*y)^2; %(norm(X(:,k),'fro'))^2 - (norm(y,'fro'))^2;
    tau = eps/(1+ eps*gn + sqrt(1 + eps*gn) );
    mu  = 1 - tau*gn;
    yp  = mu*y + tau*g;
    hp  = Z'* yp;
    e   = (tau/mu)*(Z*g - (hp'*g)*g);
    Z   = (1/beta)* (Z - g*hp' + e*g');
    ep  = mu* X(:,k) - W*yp;
    W   = W + ep*g';
    
    %% DOA tracking using ESPRIT algorithm
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
    %%
end
end