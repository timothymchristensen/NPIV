%   Computes 100(1-\alpha)% uniform confidence band for NPIV estimator via
%   bootstrap

function [xl, xu, hhat] = NPIVucb(Y,B,Psi,eps,alpha,Nboot,Psixx)

[c,~] = NPIVreg(Y,Psi,B,eps);
uhat  = Y - Psi*c;
[n,K] = size(B);
Gb    = B'*B/n;
S     = B'*Psi/n;
Buhat = B.*(uhat*ones(1,K));
Omg   = Buhat'*Buhat/n;
SGl   = inv(S'*inv(Gb)*S)*S'*inv(Gb);
Q     = SGl*Omg*SGl';
nxx   = length(Psixx);
Vxx   = zeros(nxx,1);
Sxx   = zeros(nxx,1);
for i = 1:nxx
    Vxx(i) = Psixx(i,:)*Q*Psixx(i,:)';
    Sxx(i) = sqrt(Vxx(i));
end

zalpha= NPIVucbcrit(Psixx,Sxx,SGl,B,uhat,alpha,Nboot,n);


hhat  = Psixx*c;
v     = Sxx/sqrt(n)*zalpha;
hxxal = hhat*ones(1,length(alpha));
xl    = hxxal-v;
xu    = hxxal+v;



%   bootsrap subroutine
function crit = NPIVucbcrit(Psixx,Sxx,SGl,B,uhat,alpha,Nboot,n)

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
    Rb    = (Psixx*R)./Sxx;
    zb(i) = max(abs(Rb));
end

crit = quantile(zb,alpha);