function [X,U_tr] = data_generator(n,T,r,epsilon)
% n       : dimension of observations
% T       : number of observations
% r       : rank of the underlying subspace (supposed to be fixed over time)
% epsilon : time-varying factor

% Author  : Le Trung Thanh  (University of Orleans, France)
% Email   : thanhle88.tbt@gmail.com


X = zeros(n,T);
U = randn(n,r);
for t = 1 : T
    wt        = randn(r,1);
    X(:,t)    = U * wt ;
    U_tr{1,t} = U;
    U         = U + epsilon(t)*randn(n,r);
end

end