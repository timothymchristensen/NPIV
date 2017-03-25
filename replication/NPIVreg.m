%   Estimate the NPIV model
%   Y_i = h_0(X_i) + e_i
%   inputs:     Y   = (N x 1) vector of observations on dependent variable
%               Psi = (N x J) matrix of basis functions for endog. variable
%               B   = (N x K) matrix of basis functions for exog. variable
%               eps = small ridge penalty for denominator matrix (set = 0
%                     for no penalization
%
%   outputs:    c   = coefficient vector
%               Hhat= estimator of h(x) at each data point
%               l   = minimum eigenvalue of Psi gram matrix
%               s   = minimum singular value of B'Psi/N


function [c,Hhat,l,s] = NPIVreg(Y,Psi,B,eps)

[N, C] = size(Psi);

%   sieve 2SLS
Gb = B'*B/N;
Gb1= inv(Gb);
S  = Psi'*B/N;
D  = S*Gb1*S';
d  = S*Gb1*(B'*Y/N);
c  = (D+eps*eye(C))\d;
Hhat = Psi*c;

% if nargout == 4
    Gpsi = Psi'*Psi/N;
    l    = min(eig(Gpsi));
    So   = inv(sqrtm(Gpsi))*S*inv(sqrtm(Gb));
    s    = min(svd(So));
% end