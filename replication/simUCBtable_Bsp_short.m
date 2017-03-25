%   UCB Monte Carlo Tables for CP design with B-spline basis

function simUCBtable_Bsp_short(N,Nsim,Nboot,SJ,SK,rJ,rK,design,s2)

fileID=fopen('sim_ucb_out.txt','a+');
fprintf(fileID,'%1.0f',clock);
fprintf(fileID,'\n');
fprintf(fileID,'design = %1.0f \n',design);
fprintf(fileID,'s2     = %3.4f \n',s2);
fprintf(fileID,'SJ     = %1.0f \n',SJ);
fprintf(fileID,'SK     = %1.0f \n',SK);
fprintf(fileID,'\n');

CP = simUCB_CP_Bsp_short(N,Nsim,Nboot,SJ,SK,rJ,rK,design,s2);
fprintf(fileID,'%1.0f %1.0f %3.3f %3.3f %3.3f \n',[rJ rK mean(CP)]);
            
fprintf(fileID,'\n');
fclose(fileID);