clear all
close all

fileID=fopen('Identifiable_Parm.txt');
C = textscan(fileID,'%s','delimiter','\n');
A(:,1)=strcat(C{1,1}); 
[R,C]=size(A); % R being the number of constructs 
P=[];

randn('state',sum(100*clock)); 
rand('state',sum(100*clock));

for CSTR_LV=1:R % The very outter layer of the loop, to access the textfile
    
    % Retrieve identifiable param at each construct level
    B=char(A{CSTR_LV,1}); 
    D(CSTR_LV,:)=strsplit(B,{': '});
    E=char(D{CSTR_LV,2}); 
    F{CSTR_LV,:}=strsplit(E,{' '});
    G=char(F{CSTR_LV,:}); % 
    Identifiable_Param=str2num(G)'; %This number changes every loop
   
    % Obtain DF datafile (guess) for estimation
    DF = Datafile(CSTR_LV,Identifiable_Param,P);
    % Obtain experimental datafile
    Exp_Data_DF= Exp_Data(CSTR_LV,R);
    
    clearvars P_Estimated_mtr

    [t,x] = Parameterization_Call_ODE(DF);  % solves mass balance odes @ each construct level
    GM_t = transpose(x(:,DF.ODE_size));     % stores estimation solution obtained from inital parameter guess

    for i=1:(DF.nstep+1)
        GM_t_Norm(i)=(GM_t(i)-GM_t(1))./(Exp_Data_DF.Max-Exp_Data_DF.Min); %normalize Gm_t by experimental data 
    end

    figure(CSTR_LV)   %generate a profile with initial guess, pink line
    plot(t,GM_t_Norm,'m-','LineWidth',2); 
    axis([0 7200 0 1]);
    xlabel('time'), ylabel('concentration')
    hold on

    GM_exp=Exp_Data_DF.Data;    % stores experimental data
    GM_exp_norm=GM_exp./(Exp_Data_DF.Max-Exp_Data_DF.Min);  %normalize experimental data

    plot(t,GM_exp_norm,'b--'); %Experimental lines
    hold on
    
    loopsz=DF.iter_num;
    GM_opt_Norm_sum=zeros(1,(DF.nstep+1));
    
    for i =1:(loopsz)
        i
        % Calls fmincon to find optimal P vectors to satisfy min(G_exp(t)^2-GM_t(t)^2)
        [P_Estimated,P_all,GM_opt,fval,exitflag]=Estimation_fmincon(Exp_Data_DF,DF,Identifiable_Param,i);
        P_Estimated_mtr(i,:)=P_Estimated;
        P(i,:)=P_all;
        
        for j=1:(DF.nstep+1) 
            %normalize GM_opt
            GM_opt_Norm(j)=(GM_opt(j)-GM_opt(1))./(Exp_Data_DF.Max-Exp_Data_DF.Min);       
        end
        
        GM_Mtx(i,:)=GM_opt_Norm;  
        GM_opt_Norm_sum=GM_opt_Norm_sum+GM_opt_Norm;
    end
    
    for j=1:(DF.nstep+1)
        GM_std(j)=std(GM_Mtx(:,j));  % find standard deviation of the plot
    end 

    GM_MEAN=GM_opt_Norm_sum./(loopsz); 
    GM_UB=GM_MEAN+1.96*GM_std; %95% confidence interval
    GM_LB=GM_MEAN-1.96*GM_std;


    plot(t,GM_Mtx ,'g');%individual modeled lines
    hold on
    plot(t,GM_MEAN,'k','LineWidth',2);% mean modeled line
    hold on
    plot(t,GM_UB,'k--',t,GM_LB,'k--');% 95% confidence
    hold all

    
    if CSTR_LV==R
        dlmwrite('P_solution.txt',P','delimiter','\t');  % fscanf used to analyse the data only reads through rows and put it into columns, transpose the matrix here so that it matched ths shape
    end


end

