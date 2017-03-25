function [f_CS,f_DWL,Df_CS,Df_DWL] = cs_dwl_est(p0,p1,dp,y,c,rP,sP,ktsP,Pmin,Pmax,rY,sY,ktsY,Ymin,Ymax)

pseq = (p1:-dp:p0)';
Sseq = zeros(length(pseq),1);

%   demand at p1 and y
p_init = (p1-Pmin)/(Pmax-Pmin);
y_init = (y-Ymin)/(Ymax-Ymin);
psip_init = bsplinemat(p_init,sP,rP,ktsP);
psiy_init = bsplinemat(y_init,sY,rY,ktsY);
psi_init = kron(psip_init,psiy_init);
h_init = psi_init*c;

%   solve ODE for CS
for i = 2:length(pseq)
    
    %   evaluate hhat(p(t),y-S(p(t))
    p_temp = (pseq(i)-Pmin)/(Pmax-Pmin);
    y_temp = (y-Sseq(i-1)-Ymin)/(Ymax-Ymin);
    psip_temp = bsplinemat(p_temp,sP,rP,ktsP);
    [psiy_temp,Dpsiy_temp] = bsplinemat(y_temp,sY,rY,ktsY);
    psi_temp = kron(psip_temp,psiy_temp);
    htemp = psi_temp*c;
    
    %   solve CS differential equation
    Sseq(i,1) = Sseq(i-1,1) + htemp*dp;
    
end

[pseq,ix] = sort(pseq,'ascend');

%   CS and DWL functionals
f_CS  = Sseq(ix);
f_DWL = f_CS - (p1-pseq)*h_init;

%   numerical integration step for derivatives
yseq = y-f_CS;
[psipseq,~] = bsplinemat((pseq-Pmin)/(Pmax-Pmin),sP,rP,ktsP);
[psiyseq,Dpsiyseq] = bsplinemat((yseq-Ymin)/(Ymax-Ymin),sY,rY,ktsY);

psiseq = zeros(length(pseq),(rP+sP)*(rY+sY));
Dpsiseq = zeros(length(pseq),(rP+sP)*(rY+sY));

for i = 1:length(pseq)
    psiseq(i,:) = kron(psipseq(i,:),psiyseq(i,:));
    Dpsiseq(i,:) = kron(psipseq(i,:),Dpsiyseq(i,:))/(Ymax-Ymin);
end

hseq = psiseq*c;
Dhseq = Dpsiseq*c;

Dfseq = zeros(length(pseq),(rP+sP)*(rY+sY));

for i = 1:length(pseq)
    
    Dtemp = zeros(1,(rP+sP)*(rY+sY));
    
    for j = i:length(pseq)
        
        intexp = exp(-sum(Dhseq(i:j)*dp));
        
        Dtemp = Dtemp + psiseq(j,:)*intexp*dp;
        
    end
    
    Dfseq(i,:) = Dtemp;
    
end

%   CS and DWL derivative terms
Df_CS = Dfseq;
Df_DWL = Df_CS - (p1-pseq)*psi_init;



