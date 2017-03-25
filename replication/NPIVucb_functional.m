%   Computes 100(1-\alpha)% uniform confidence band for NPIV estimator via
%   bootstrap

function [xl, xu] = NPIVucb_functional(Y,B,Psi,eps,alpha,Nboot,f,Df)

[c,~] = NPIVreg(Y,Psi,B,eps);
uhat  = Y - Psi*c;
[n,K] = size(B);
Gb    = B'*B/n;
S     = B'*Psi/n;
Buhat = B.*(uhat*ones(1,K));
Omg   = Buhat'*Buhat/n;
SGl   = inv(S'*inv(Gb)*S)*S'*inv(Gb);
Q     = SGl*Omg*SGl';
nxx   = length(f);
Vxx   = zeros(nxx,1);
Sxx   = zeros(nxx,1);
for i = 1:nxx
    Vxx(i) = Df(i,:)*Q*Df(i,:)';
    Sxx(i) = sqrt(Vxx(i));
end

zalpha= NPIVucbcrit_functional(Df,Sxx,SGl,B,uhat,alpha,Nboot,n);

v	= Sxx/sqrt(n)*zalpha;
fal = f*ones(1,length(alpha));
xl	= fal-v;
xu	= fal+v;



%   bootsrap subroutine
function crit = NPIVucbcrit_functional(Df,Sxx,SGl,B,uhat,alpha,Nboot,n)

zb   = zeros(Nboot,1);

for i = 1:Nboot
    %  N(0,1)
%     w     = randn(n,1);
    %  Re-centered exponential
%     w     = exprnd(1,n,1)-1;
    % Rademacher
%     w     = 2*(rand(n,1) > 0.5)-1;
    % Mammen two-point
    w     = 0.5 - 0.5*sqrt(5)*(2*(rand(n,1) < (sqrt(5)+1)/(2*sqrt(5)))-1);
    R     = SGl*(B'*(uhat.*w))/sqrt(n);
    Rb    = (Df*R)./Sxx;
    zb(i) = max(abs(Rb));
end

crit = quantile(zb,alpha);