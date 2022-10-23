function [theta] = DOA_FAPI(X,beta,r)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% r     : target rank (number of signals)

[n, N] = size(X);
W      = eye(n,r);
Z      = eye(r);

theta = [];

for k = 1:N
    
    y  = W'*X(:,k);
    h  = Z*y;
    g  = h/(beta + y'*h);
    gn = (norm(g,'fro'))^2;

    eps = (norm(X(:,k),'fro'))^2 - (norm(y,'fro'))^2;
    tau = eps/(1+ eps*gn + sqrt(1 + eps*gn) );
    mu  = 1 - tau*gn;
    yp  = mu*y + tau*g;
    hp  = Z'* yp;
    e   = (tau/mu)*(Z*g - (hp'*g)*g);
    Z   = (1/beta)* (Z - g*hp' + e*g');
    ep  = mu* X(:,k) - W*yp;
    W   = W + ep*g';
    
    % DOA tracking using ESPRIT algorithm
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