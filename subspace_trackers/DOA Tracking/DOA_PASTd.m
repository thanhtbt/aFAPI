function [theta,rho,eta,eigvalue] = DOA_PASTd(X,bet,r)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% r     : target rank (number of signals)


[n, N] = size(X);
%Initialization
W = eye(n,r);
d = ones(r,1);
V = orth(X);
rho = zeros(1,N);
eta = zeros(1,N);
theta = [];
%Processing
for k = 1:N
    
    xt = X(:,k);
    for j = 1 : r
        y = (W(:,j))' * xt;
        d(j) = bet * d(j) + abs(y)^2;
        e = xt - W(:,j)*y;
        W(:,j) = W(:,j) + e * (y'/d(j));
        xt = xt - W(:,j)*y;
    end
       
%   % performance of subspace tracking
    rho(k)=  abs(trace(W'*(eye(n)-V*V')*W)/trace(W'*(V*V')*W));
    eta(k)=  norm(W' * W - eye(r))^2;
    
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
    eigvalue = d;
end