function [theta] = DOA_TRPAST(X,beta,alpha,r)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% alpha : divergence

% Author: Le Trung Thanh  (University of Orleans, France)
% Email : thanhle88.tbt@gmail.com
% Cite  : L.T. Thanh, A.M. Rekavandi, S. Abd-Krim, & K. Abed-Meraim.
          ... “Robust Subspace Tracking With Contamination via Alpha-Divergence”. 
          ... ICASSP 2023 (submitted).


[n, N] = size(X);
theta = [];

W   = eye(n,r);
P   = eye(r);
c   = (1-alpha)/2;

Pt = P;
%Processing
for k = 1 : N
    %% TRPAST section
    xt    = X(:,k);
    % Calculate weight
    yt    = W'*xt;
    et    = xt - W*yt;
    del_t = norm(et)^2;
    wt    = exp(-c*del_t);  
    % Update subspace
    ht    = Pt * yt;
    gt    = wt*ht/(beta + wt*yt'*ht);
    Pt    = beta^(-1) * ( (Pt - gt*ht') );
    W     = W + et * gt';
        
    
    %% DOA tracking using ESPRIT algorithm
    D1 = W(1:n-1,:);
    D2 = W(2:n,:);
    
    % TLS solution
    [~, ~, U] = svd([D1 D2]);
    U12 = U(1:r,(r+1):(2*r));
    U22 = U((r+1):(2*r),(r+1):(2*r));
    Phi = -U12*pinv(U22); % TLS solution for Psi
    
    [~,Dest] = eigs(Phi);  
    west     = angle((diag(Dest)))/(2*pi);
    
    theta = [theta sort(west)];

end