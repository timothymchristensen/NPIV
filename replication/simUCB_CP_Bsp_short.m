function CP = simUCB_CP_Bsp_short(N,Nsim,Nboot,SJ,SK,rJ,rK,design,s2)

%   Preliminaries
rng(7654321);
dx   = 1e-4;
x    = (.05:dx:.95)';
eps  = 1e-5;
alpha= [0.90 0.95 0.99];
nal  = length(alpha);

%   Function
H0   = H0fun(x,design);

%   Generate bases
PSI = bsplinemat(x,SJ,rJ);

%   Preallocate
CP   = zeros(Nsim,nal);

%   Loop
for n = 1:Nsim
    %   Simulate
    [Y,X,W,~] = simNPIV(N,s2,design); 
    
    Psi  = bsplinemat(X,SJ,rJ);
    B    = bsplinemat(W,SK,rK);
    
    %   Uniform confidence band coverage
    [xl, xu, ~] = NPIVucb(Y,B,Psi,eps,alpha,Nboot,PSI);
    for i = 1:nal
        CP(n,i) = min( xl(:,i) < H0)*min( H0 < xu(:,i) );
    end
    
end




