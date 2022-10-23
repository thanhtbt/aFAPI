function [Ut,rho, eta] = TRPAST(X,beta,alpha,U_tr)
% X    : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta : forgetting factor
% alpha: divergence
% U_tr : the set of true subspaces with time

% Author: Le Trung Thanh  (University of Orleans, France)
% Email : thanhle88.tbt@gmail.com
% Cite  : L.T. Thanh, A.M. Rekavandi, S. Abd-Krim, & K. Abed-Meraim.
          ... “Robust Subspace Tracking With Contamination via Alpha-Divergence”.
          ... ICASSP 2023 (submitted).
% Source paper: A.M. Rekavandi, S. Abd-Krim, & K. Abed-Meraim, 
               ... “TRPAST: A tunable and robust projection approximation subspace tracking method,” 
               ... IEEE Trans. Signal Process., 2022 (submitted).
[n, N] = size(X);
r   = size(U_tr{1,1},2);
W   = eye(n,r);
P   = eye(r);
c   = (1-alpha)/2;


rho = zeros(1,N);
eta = zeros(1,N);

Pt = P;
%Processing
for k = 1 : N
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
    
    
    %% Performance Evaluation
    
    % Performance Evaluation
    V       = orth(U_tr{1,k}); % true subspace at time t
    rho(k)  = abs(trace(pinv(W)*(eye(n)-V*V')*W)/trace(pinv(W)*(V*V')*W));
    eta(k)  = sin(subspace(W,V));
    Ut{1,k} = W;
end