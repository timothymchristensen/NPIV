% Master file for running UCB Engel curve empirical application

clear all, clc

% Load the data
data = xlsread('eng95a.xls');
sfood = data(:,3);    % food in
scater = data(:,4);   % food out
salc = data(:,5);     % alcohol
sfuel = data(:,6);    % fuel
smotor = data(:,7);   % travel
sfares = data(:,8);   % fares
slserv = data(:,9);   % leisure goods
lnx = data(:,10);     % log nondurable expenditure
lngw = data(:,12);    % log gross wage (instrument)
kids = data(:,14);    % kids = 1

Ymat = [sfood,scater,salc,sfuel,smotor,slserv];
names = ['food in ';...
         'food out';...
         'alcohol ';...
         'fuel    ';...
         'travel  ';...
         'leisure '];

%   dimension of sieve basis
rJ   = 5;
rK   = 5;
SJ   = 0;
SK   = 4;

%   Preliminaries
rng(7654321);
dx   = 1e-5;
eps  = 1e-5;
alpha= [0.90 0.95];
nal  = length(alpha);
Smax = length(SJ);
Nboot= 1e3;


%   engel curves -- kids
ix = logical(kids);
lnxk  = lnx(ix);
lngwk = lngw(ix);
X = (lnxk-min(lnxk))/(max(lnxk)-min(lnxk));
W = (lngwk-min(lngwk))/(max(lngwk)-min(lngwk));
x = (0:dx:1)';
x0= 4.75+x*(6.25-4.75);
x1= (x0-min(lnxk))/(max(lnxk)-min(lnxk));
for i = 1:6
    
    Y   = Ymat(ix,i);
    
    [Psi,~,~,ktsJ] = bsplinemat_quantile(X,SJ,rJ);
    [PSI, DPSI] = bsplinemat(x1,SJ,rJ,ktsJ);
    
    B   = bsplinemat_quantile(W,SK,rK);
    
    [xl, xu, hhat] = NPIVucb(Y,B,Psi,eps,alpha,Nboot,PSI);
    
    figure(1)
    subplot(2,3,i)
    plot(x0,hhat,'r','LineWidth',1.2)
    hold on
    plot(x0,xl(:,1),'k-','LineWidth',1.5)
    plot(x0,xu(:,1),'k-','LineWidth',1.5)
    plot(x0,xl(:,2),'k--','LineWidth',1.5)
    plot(x0,xu(:,2),'k--','LineWidth',1.5)
    axis([4.7 6.3 0.0 0.4])
    ylabel(names(i,:))
    grid on
    hold off
    
    [xl, xu, hhat] = NPIVucb(Y,B,Psi,eps,alpha,Nboot,DPSI);
    
    figure(2)
    subplot(2,3,i)
    plot(x0,hhat,'r','LineWidth',1.2)
    hold on
    plot(x0,xl(:,1),'k-','LineWidth',1.5)
    plot(x0,xu(:,1),'k-','LineWidth',1.5)
    plot(x0,xl(:,2),'k--','LineWidth',1.5)
    plot(x0,xu(:,2),'k--','LineWidth',1.5)
    axis([4.7 6.3 -2 3])
    ylabel(names(i,:))
    grid on
    hold off
    
end


