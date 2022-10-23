function [Ut,rho, eta] = PASTd_Tracking(X,beta,U_tr)
% X    : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta : forgetting factor
% U_tr : the set of true subspaces with time

%Initialization
[n, N] = size(X);
r      = min(size(U_tr{1,1}));
W      = eye(n,r);
Z      = eye(r);
d      = ones(r,1);


rho = zeros(1,N);
eta = zeros(1,N);

%Processing
for k = 1:N
    
    xt = X(:,k);
    for j = 1 : r
        y      = (W(:,j))' * xt;
        d(j)   = beta * d(j) + abs(y)^2;
        e      = xt - W(:,j)*y;
        W(:,j) = W(:,j) + e * (y'/d(j));
        xt     = xt - W(:,j)*y;
    end
    % Performance Evaluation
    V       = orth(U_tr{1,k}); % true subspace at time t
    rho(k)  = abs(trace(pinv(W)*(eye(n)-V*V')*W)/trace(pinv(W)*(V*V')*W));
    eta(k)  = sin(subspace(W,V));
    Ut{1,k} = W;
end