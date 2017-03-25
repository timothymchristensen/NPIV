% Master file for running UCB CS DWL empirical application

clear all, clc

% Load the data
data = csvread('bhp_new.csv',2,1);
Q = data(:,1);      % quantity
P = data(:,2);      % price
Y = data(:,3);      % income
W = data(:,4);      % distance instrument
hhsize = data(:,5);	% household size
driver = data(:,6);	% number of drivers

%   dimension of sieve basis
rP   = 5;
rY   = 5;
rW   = 5;
sP   = 0;
sY   = 0;
sW   = 3;

%   Preliminaries
rng(7654321);
dx   = 1e-5;
eps  = 1e-5;
alpha= [0.90 0.95];
nal  = length(alpha);
Nboot= 1e3;

%   scale regressors and instruments
ix = (Y>25000&Y<100000&hhsize<=6&driver<=2&Q<=10000);
Q = Q(ix);
P = P(ix);
Y = Y(ix);
W = W(ix);
p = (P-min(P))/(max(P)-min(P));
y = (Y-min(Y))/(max(Y)-min(Y));
w = (W-min(W))/(max(W)-min(W));
x = (0:dx:1)';

% univariate basis functions
[PsiP,~,~,ktsP] = bsplinemat_quantile(p,sP,rP);
[PsiY,~,~,ktsY] = bsplinemat_quantile(y,sY,rY);
[PsiW,~,~,ktsW] = bsplinemat_quantile(w,sW,rW);

% tensor product
Psi = zeros(length(p),(rP+sP)*(rY+sY));
B   = zeros(length(p),(rW+sW)*(rY+sY));
for i = 1:length(p)
    Psi(i,:) = kron(PsiP(i,:),PsiY(i,:));
    B(i,:)= kron(PsiW(i,:),PsiY(i,:));
end

% estimate demand function by NPIV
[c,~] = NPIVreg(Q,Psi,B,eps);

% estimate demand function by LS
c_ols = inv(Psi'*Psi)*Psi'*Q;

% CS and DWL estimation
y = 72500;
p0= 1.20;
p1= 1.40;
dp= 1e-4;
Pmin = min(P);
Pmax = max(P);
Ymin = min(Y);
Ymax = max(Y);
[f_CS,f_DWL,Df_CS,Df_DWL] = cs_dwl_est(p0,p1,dp,y,c,rP,sP,ktsP,Pmin,Pmax,rY,sY,ktsY,Ymin,Ymax);

f_CS_NPIV_high = f_CS;
f_DWL_NPIV_high = f_DWL;

% UCB for CS and DWL 
[xl_CS, xu_CS]   = NPIVucb_functional(Q,B,Psi,eps,alpha,Nboot,f_CS,Df_CS);
[xl_DWL, xu_DWL] = NPIVucb_functional(Q,B,Psi,eps,alpha,Nboot,f_DWL,Df_DWL);

figure(1)
subplot(2,2,1)
pseq = (p0:dp:p1)';
plot(pseq,f_CS,'r','LineWidth',1.2)
hold on
plot(pseq,xl_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_CS(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_CS(:,2),'k--','LineWidth',1.5)
axis([p0 p1 0 350])
xlabel('initial price ($/gal)')
ylabel('CS ($)')
grid on
hold off
subplot(2,2,3)
pseq = (p0:dp:p1)';
plot(pseq,f_DWL,'r','LineWidth',1.2)
hold on
plot(pseq,xl_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_DWL(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_DWL(:,2),'k--','LineWidth',1.5)
axis([p0 p1 -50 120])
xlabel('initial price ($/gal)')
ylabel('DL ($)')
grid on
hold off

% repeat for lower income level
y = 42500;
[f_CS,f_DWL,Df_CS,Df_DWL] = cs_dwl_est(p0,p1,dp,y,c,rP,sP,ktsP,Pmin,Pmax,rY,sY,ktsY,Ymin,Ymax);

f_CS_NPIV_low = f_CS;
f_DWL_NPIV_low = f_DWL;

% UCB for CS and DWL 
[xl_CS, xu_CS]   = NPIVucb_functional(Q,B,Psi,eps,alpha,Nboot,f_CS,Df_CS);
[xl_DWL, xu_DWL] = NPIVucb_functional(Q,B,Psi,eps,alpha,Nboot,f_DWL,Df_DWL);

figure(1)
subplot(2,2,2)
pseq = (p0:dp:p1)';
plot(pseq,f_CS,'r','LineWidth',1.2)
hold on
plot(pseq,xl_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_CS(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_CS(:,2),'k--','LineWidth',1.5)
axis([p0 p1 0 350])
xlabel('initial price ($/gal)')
ylabel('CS ($)')
grid on
hold off
subplot(2,2,4)
pseq = (p0:dp:p1)';
plot(pseq,f_DWL,'r','LineWidth',1.2)
hold on
plot(pseq,xl_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_DWL(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_DWL(:,2),'k--','LineWidth',1.5)
axis([p0 p1 -50 120])
xlabel('initial price ($/gal)')
ylabel('DL ($)')
grid on
hold off

fprintf('mean & %3.0f & %3.2f & %3.0f \\\\ \n',mean([Q P Y]))
fprintf('25th percentile & %3.0f & %3.2f & %3.0f \\\\ \n',quantile([Q P Y],0.25))
fprintf('median & %3.0f & %3.2f & %3.0f \\\\ \n',median([Q P Y]))
fprintf('75th percentile & %3.0f & %3.2f & %3.0f \\\\ \n',quantile([Q P Y],0.75))
fprintf('std dev & %3.0f & %3.2f & %3.0f \\\\ \n',std([Q P Y]))


% repeat for series LS higher income level
y = 72500;
[f_CS,f_DWL,Df_CS,Df_DWL] = cs_dwl_est(p0,p1,dp,y,c_ols,rP,sP,ktsP,Pmin,Pmax,rY,sY,ktsY,Ymin,Ymax);

% UCB for CS and DWL 
[xl_CS, xu_CS]   = NPIVucb_functional(Q,Psi,Psi,eps,alpha,Nboot,f_CS,Df_CS);
[xl_DWL, xu_DWL] = NPIVucb_functional(Q,Psi,Psi,eps,alpha,Nboot,f_DWL,Df_DWL);

figure(2)
subplot(2,2,1)
pseq = (p0:dp:p1)';
plot(pseq,f_CS,'r','LineWidth',1.2)
hold on
plot(pseq,xl_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_CS(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_CS(:,2),'k--','LineWidth',1.5)
plot(pseq,f_CS_NPIV_high,'r-.','LineWidth',1.2)
axis([p0 p1 0 350])
xlabel('initial price ($/gal)')
ylabel('CS ($)')
grid on
hold off
subplot(2,2,3)
pseq = (p0:dp:p1)';
plot(pseq,f_DWL,'r','LineWidth',1.2)
hold on
plot(pseq,xl_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_DWL(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_DWL(:,2),'k--','LineWidth',1.5)
plot(pseq,f_DWL_NPIV_high,'r-.','LineWidth',1.2)
axis([p0 p1 -50 120])
xlabel('initial price ($/gal)')
ylabel('DL ($)')
grid on
hold off

% repeat for series LS lower income level
y = 42500;
[f_CS,f_DWL,Df_CS,Df_DWL] = cs_dwl_est(p0,p1,dp,y,c_ols,rP,sP,ktsP,Pmin,Pmax,rY,sY,ktsY,Ymin,Ymax);

% UCB for CS and DWL 
[xl_CS, xu_CS]   = NPIVucb_functional(Q,Psi,Psi,eps,alpha,Nboot,f_CS,Df_CS);
[xl_DWL, xu_DWL] = NPIVucb_functional(Q,Psi,Psi,eps,alpha,Nboot,f_DWL,Df_DWL);

figure(2)
subplot(2,2,2)
pseq = (p0:dp:p1)';
plot(pseq,f_CS,'r','LineWidth',1.2)
hold on
plot(pseq,xl_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_CS(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_CS(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_CS(:,2),'k--','LineWidth',1.5)
plot(pseq,f_CS_NPIV_low,'r-.','LineWidth',1.2)
axis([p0 p1 0 350])
xlabel('initial price ($/gal)')
ylabel('CS ($)')
grid on
hold off
subplot(2,2,4)
pseq = (p0:dp:p1)';
plot(pseq,f_DWL,'r','LineWidth',1.2)
hold on
plot(pseq,xl_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xu_DWL(:,1),'k-','LineWidth',1.5)
plot(pseq,xl_DWL(:,2),'k--','LineWidth',1.5)
plot(pseq,xu_DWL(:,2),'k--','LineWidth',1.5)
plot(pseq,f_DWL_NPIV_low,'r-.','LineWidth',1.2)
axis([p0 p1 -50 120])
xlabel('initial price ($/gal)')
ylabel('DL ($)')
grid on
hold off
