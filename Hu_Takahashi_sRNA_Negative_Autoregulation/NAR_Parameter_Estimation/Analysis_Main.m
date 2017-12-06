clear all
close all
CSTR_LV=1;
DF = Datafile_Analysis(CSTR_LV);
R=DF.Total_constructs;

for CSTR_LV=1:R  

    DF = Datafile_Analysis(CSTR_LV);    %Access solved estimation datafile from P_solutions
    Exp_Data_DF= Exp_Data(CSTR_LV,R);   %Access experimental datafile

    GM_exp=Exp_Data_DF.Data;
    GM_exp_norm=GM_exp./(Exp_Data_DF.Max-Exp_Data_DF.Min);

    loopsz=DF.iter_num;
    GM_t_Norm_sum=zeros(1,(DF.nstep+1));
    
    for i =1:(loopsz)
        i
        P_lib=DF.Parameter_library;
        NParameters=DF.Num_Parameters;
% -----------------------------------------   % drawing across columns      
        for j=1:(NParameters)        
          index=randi(loopsz,1);
          P(j)=P_lib(j,index);
        end
% ------------------------------------------  end        
% %------------------------------------------  % drawing by columns
%         index=randi((loopsz),1);  
%         P=P_lib(index,:);
%------------------------------------------  end   
% % ------------------------------------------  % drawing by order
%           P=P_lib(i,:);
% ------------------------------------------  end   
        DF.Initial_Parameters=P;
        [t,x] = Parameterization_Call_ODE(DF);%solves mass balance odes @ each construct level
        GM_t= transpose(x(:,DF.ODE_size));
        
        for j=1:(DF.nstep+1)
            GM_t_Norm(j)=(GM_t(j)-GM_t(1))./(Exp_Data_DF.Max-Exp_Data_DF.Min); 
        end
    
        GM_Mtx(i,:)=GM_t_Norm;  
        GM_t_Norm_sum=GM_t_Norm_sum+GM_t_Norm;    
    end

    for j=1:(DF.nstep+1)
        GM_std(j)=std(GM_Mtx(:,j));
    end 
% 

    GM_MEAN=GM_t_Norm_sum./(loopsz); 
    GM_UB=GM_MEAN+1.96*GM_std;
    GM_LB=GM_MEAN-1.96*GM_std;

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
    
    DATA=[GM_MEAN;GM_UB;GM_LB;GM_exp_norm']; % Saved text files contain mean,lower bound, upper bound and normalized experimenl data.
    filename=['PLOT_DATA_CSTR_',num2str(CSTR_LV),'.txt'];
    dlmwrite(filename,DATA','delimiter','\t');

end

