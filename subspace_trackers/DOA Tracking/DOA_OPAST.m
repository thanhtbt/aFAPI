function [theta,rho,eta] = DOA_OPAST(X,beta,r)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% r     : target rank (number of signals)


[n, N] = size(X);
%Initialization
V = orth(X);
W   = eye(n,r);
Z   = eye(r);
% rho = zeros(1,N);
% eta = zeros(1,N);
theta = [];
%Processing
for k = 1:N
    y     = W'*X(:,k);
    q     = Z*y/beta;
    GAMA  = 1/(1 + y'*q);
    TAUX  = 1/(norm(q,'fro'))^2*(1/sqrt(1 + (norm(q,'fro'))^2*GAMA^2*(norm(X(:,k),'fro')^2-norm(y,'fro')^2 )) - 1);
    ee    = W*(TAUX*q - GAMA*(1 + TAUX*(norm(q,'fro'))^2)*y) + (1 + TAUX*(norm(q,'fro'))^2)*GAMA*X(:,k);
    Z     = Z/beta - GAMA*(q*q');
    W     = W + ee*q';
    
%     % performance of subspace tracking
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

end