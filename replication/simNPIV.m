function [Y,X,W,H0] = simNPIV(N,s2,design)

Sig  = [1 .5 0; .5 1 0; 0 0 1];
R    = chol(Sig);
Z    = randn(N,3)*R;

U    = Z(:,1);
V    = Z(:,2);
W    = Z(:,3);
X    = normcdf((W + V)/sqrt(2));
W    = normcdf(Z(:,3));

H0   = H0fun(X,design);
Y    = H0 + sqrt(s2)*U;

