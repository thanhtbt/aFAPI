function [theta] = DOA_LORAF(X,beta,r)
% X     : an (n x N) data matrix collecting the N observation (n x 1) vectors
% beta  : forgetting factor (0 < beta<= 1)
% r     : target rank (number of signals)

[n, N] = size(X);
theta = [];

%Initialization
D   = eye(n,r);

S = 0.00001*eye(r);
v = zeros(r,1);

% performance assessment
rho = zeros(1,N);
eta = zeros(1,N);
%Processing

for k = 1:N
    z = X(:,k);
    h = D'*z;
    Z = z'*z - h'*h;
    if Z == 0
        print('Oh la la');
    end
    u = S*v;
    XX = beta*S - 2*beta*u*v' + h*h';
    b = sqrt(Z)*pinv(XX')*h;
    bet = 4*(b'*b + 1);
    p2 = 1/2 + 1/sqrt(bet);
    gam = (1-2*p2)/(2*sqrt(p2));
    del = sqrt(p2)/sqrt(Z);
    v = gam*b;
    S = XX-(1/del)*v*h';
    w = del*h-v;
    e = del*z - D*w;
    D = D - 2*e*v';
    
    % DOA tracking using ESPRIT algorithm
    D1 = D(1:n-1,:);
    D2 = D(2:n,:);
    
    % TLS solution
    [~, ~, U] = svd([D1 D2]);
    U12 = U(1:r,(r+1):(2*r));
    U22 = U((r+1):(2*r),(r+1):(2*r));
    Phi = -U12*pinv(U22); % TLS solution for Psi
    
    [~,Dest] = eigs(Phi);  
    west = angle((diag(Dest)))/(2*pi);
    
    theta = [theta sort(west)];
end

