function [P_Estimated,P_all,GM_opt,fval,exitflag] = Estimation_fmincon(Exp_Data_DF,DF,Identifiable_Param,i)
%This function uses its nested fmincon function to estimate identifiable
%parameters using experimental data

    CSTR_LV=DF.Construct;
    NParameters=DF.Num_Parameters;
    iter_num=DF.iter_num;
    tspan=DF.tspan;
    ODE_size=DF.ODE_size;
    P=DF.Parameter_library(i,:); %Columns are the different parameters, rows are different sets
        
    function [Error] = fminconfunc(y)
             
        timesize=Exp_Data_DF.timestep; % get experimental data 
        
        for k=1:length(Identifiable_Param)  %update P with identifiable parameters
            P(Identifiable_Param(k))=y(k);
        end
            
        [t,x] = Parameterization_Call_ODE(DF);
        GM_est=transpose(x(:,ODE_size)); %Basing estimate error on GFP prediction
        Error=0; 
           
         for j=1:(timesize)
             Error=Error+(GM_est(j)-Exp_Data_DF.avg(j)).^2; %Sum of Error 
         end
   
    end

    
    y0=P(Identifiable_Param(1));
%     [ A,b ] = fmincon_constrains(pset,P);
    for m=2:length(Identifiable_Param)
        y0=[y0,P(Identifiable_Param(m))];
    end

    lb=0.33*y0; 
    ub=3.33*y0;

    
%     [y,fval,exitflag] = fmincon(@fminconfunc,y0,A,b,[],[],lb,ub);
    [y,fval,exitflag] = fmincon(@fminconfunc,y0,[],[],[],[],lb,ub);
    P_Estimated=y;

    P_all=P;
    for i=m:(length(Identifiable_Param))
        P_all((Identifiable_Param(m)))=P_Estimated(m);
    end
    
    P=P_all;
    DF.Initial_Parameters=P;
    [t,x] = Parameterization_Call_ODE(DF);
    GM_opt=transpose(x(:,ODE_size));    
    
end


