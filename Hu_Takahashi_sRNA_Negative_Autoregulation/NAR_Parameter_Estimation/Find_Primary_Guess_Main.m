clear all
close all

CSTR_LV=1;
DF = Datafile_Analysis(CSTR_LV); %define struct
iter_num=DF.iter_num;
R=DF.Total_constructs;
Sum_error=zeros(iter_num,1);

for CSTR_LV=1:R  

    DF = Datafile_Analysis(CSTR_LV);    %Access solved estimation datafile from P_solutions
    Exp_Data_DF= Exp_Data(CSTR_LV,R);   %Access experimental datafile

    GM_exp=Exp_Data_DF.Data; 
    GM_exp_norm=GM_exp./(Exp_Data_DF.Max-Exp_Data_DF.Min); %normalize experimental data 

    GM_t_Norm_sum=zeros(1,(DF.nstep+1));  
    Norm_Exp_Avg=Exp_Data_DF.avg./Exp_Data_DF.Max; % the avarage of each time piont
    GM_error=zeros(iter_num,1);
    
    for i =1:iter_num
        i
        P_lib=DF.Parameter_library;
        NParameters=DF.Num_Parameters;
       
          P=P_lib(:,i);    % drawing by order

        DF.Initial_Parameters=P;   
        [t,x] = Parameterization_Call_ODE(DF);   %solves mass balance odes @ each construct level
        GM_t= transpose(x(:,DF.ODE_size));
        
        for j=1:(DF.nstep+1)
            GM_t_Norm(j)=(GM_t(j)-GM_t(1))./(Exp_Data_DF.Max-Exp_Data_DF.Min); %normalized solution
            GM_error(i)= GM_error(i)+abs(GM_t_Norm(j)- Norm_Exp_Avg(j)); 
        end
    
        GM_Mtx(i,:)=GM_t_Norm;  
        GM_t_Norm_sum=GM_t_Norm_sum+GM_t_Norm;    % save a sum of the normalized solution
    end

    for j=1:(DF.nstep+1)
        GM_std(j)=std(GM_Mtx(:,j));  
    end 
% 

    GM_MEAN=GM_t_Norm_sum./(iter_num); 
    GM_UB=GM_MEAN+1.96*GM_std; % 95% confidence interval would be plotted
    GM_LB=GM_MEAN-1.96*GM_std;
    
    Sum_error=Sum_error+GM_error;   %save the sum of errors 

    figure(CSTR_LV)   %generate a profile 
    plot(t,GM_Mtx ,'g');
   
    hold on
    plot(t,GM_MEAN,'k','LineWidth',3);
    plot(t,GM_UB,'k--',t,GM_LB,'k--');
    axis([0 7200 0 1]);
    xlabel('time'), ylabel('concentration')
    hold all
    
    plot(t,GM_exp_norm,'b--'); 
    hold on
    
end

[C,Index]=min(Sum_error); %Pick the guess set that generated a simulation that's closest to the experimental one
Best_fit=P_lib(:,Index)';
dlmwrite('guess.txt',Best_fit','delimiter','\t');  %generate a new guess based on previous trail for the next round to start with.