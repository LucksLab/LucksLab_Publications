clear all
close all
CSTR_LV=1;
DF = Datafile_Prediction(CSTR_LV);
R=DF.Total_constructs;

for CSTR_LV=1:R % 

    [DF] = Datafile_Prediction(CSTR_LV);
    Exp_Data_DF= Exp_Data_NAR(CSTR_LV,R);
    loopsz=DF.iter_num;
    GM_t_sum=zeros(1,(DF.nstep+1));
     for i =1:loopsz
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
%         P=P_lib(:,index);
%------------------------------------------  end   
% % ------------------------------------------  % drawing by order
%           P=P_lib(i,:);
% ------------------------------------------  end   
       
        DF.Initial_Parameters=P;
        [t,x] = NAR_Call_ODE(DF); %solves mass balance odes @ each construct level
        GM_t= transpose(x(:,DF.ODE_size));
       
    
        GM_Mtx(i,:)=GM_t;  
        GM_t_sum=GM_t_sum+GM_t;    
    end

    for j=1:(DF.nstep+1)
        GM_std(j)=std(GM_Mtx(:,j));
    end 
% 

    GM_MEAN=GM_t_sum./(loopsz); 
    GM_UB=GM_MEAN+1.96*GM_std;
    GM_LB=GM_MEAN-1.96*GM_std;

%     figure (CSTR_LV) %generate a profile  `
%     plot(t,GM_Mtx ,'g');
%    

    if CSTR_LV==1
    figure (1)
    hold all
    plot(t,GM_MEAN,'k','LineWidth',3);
    plot(t,GM_UB,'k--',t,GM_LB,'k--');
    plot(t,Exp_Data_DF.Data,'c');
    elseif CSTR_LV==2
    figure (1)
    plot(t,GM_MEAN,'b','LineWidth',3);
    plot(t,GM_UB,'b--',t,GM_LB,'b--');
    plot(t,Exp_Data_DF.Data,'r');
    elseif CSTR_LV==3
    figure (2)
    hold all
    plot(t,GM_MEAN,'g','LineWidth',3);
    plot(t,GM_UB,'g--',t,GM_LB,'g--');
    plot(t,Exp_Data_DF.Data,'c');
    elseif CSTR_LV==4
    figure (2)
    plot(t,GM_MEAN,'m','LineWidth',3);
    plot(t,GM_UB,'m--',t,GM_LB,'m--');
    plot(t,Exp_Data_DF.Data,'r');
    end
    
    xlabel('time'), ylabel('concentration')
    hold all

    
	filename=['Thalfplot_simulated_',num2str(CSTR_LV),'.txt'];
	dlmwrite(filename,GM_Mtx','delimiter','\t');

end

