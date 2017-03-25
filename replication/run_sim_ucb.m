%   Master file for running simulation experiments
%	Prints output to sim_ucb_out.txt

N    = 1e3;
Nsim = 1e3;
Nboot= 1e3;

%   Design 1

s2   = 1;
design=1;

simUCBtable_Bsp_short(N,Nsim,Nboot,1,1,4,4,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,1,2,4,4,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,1,0,4,5,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,1,1,4,5,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,0,0,5,5,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,0,1,5,5,design,s2)

%   Design 2, SJ = SK

design=2;

simUCBtable_Bsp_short(N,Nsim,Nboot,1,1,4,4,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,1,2,4,4,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,1,0,4,5,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,1,1,4,5,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,0,0,5,5,design,s2)
simUCBtable_Bsp_short(N,Nsim,Nboot,0,1,5,5,design,s2)


