function [Y,X,W,H0] = simNPIV(N,s2,design)

% Sig  = [1 .8 .5; .8 1 0; .5 0 1];
% R    = chol(Sig);
% Z    = randn(N,3)*R;
% 
% X    = normcdf(Z(:,1));
% W    = normcdf(Z(:,2));
% U    = sqrt(s2)*Z(:,3);
% H0   = H0fun(X,design);
% Y    = H0 + U;


Sig  = [1 .5 0; .5 1 0; 0 0 1];
R    = chol(Sig);
Z    = randn(N,3)*R;

U    = Z(:,1);
V    = Z(:,2);
W    = Z(:,3);
X    = normcdf((W + V)/sqrt(2));
W    = normcdf(Z(:,3));

% 
% 
% X    = normcdf(Z(:,1));
% W    = normcdf(Z(:,2));
% U    = sqrt(s2)*Z(:,3);
H0   = H0fun(X,design);
Y    = H0 + sqrt(s2)*U;

