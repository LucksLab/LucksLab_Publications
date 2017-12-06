
%This function uses numerical method to find SSM and normalize it. 
%Both the Jacobian Matrices and the P Matrices are estimated using 4th
%order central difference

function [ SSM ] = Parameterization_Numerical_SSM( DF,x,IDP)

    P=DF.Initial_Parameters;  
    t_inc=DF.t_inc; 
    nstep=DF.nstep;
    P_size=DF.Num_Parameters; 
    ODE_size=DF.ODE_size;
    CSTR_LV=DF.Construct
    z0=zeros(ODE_size,1);  % Initial value for z thats used in ode15s 

    SSM=zeros((nstep),(P_size)); 

    for k=2:(nstep+1) %Solves SSM at each time point (total of 21)
        print=['solving for SSM at point ',num2str(k),'/',num2str(nstep+1)];
        disp(print);
        for j=1:(P_size) %Solves dx/dp for all x and all p at time point k
            tspan=0:(t_inc):((k-1)*t_inc);
            [t,x] = Parameterization_Call_ODE(DF); %solve for all Xs in timeframe set by tspan
        
            [ J_mtrx ] = Parameterization_Jacobian( P,x,k,CSTR_LV,ODE_size ); % get the jacobian matrix
        
            [ P_mtrx ] = Parameterization_Pmatrix (P,x,k,CSTR_LV,ODE_size,P_size,j);  % get the pmatrix
        
            [t,z]=ode15s(@senfunction,tspan,z0,[],J_mtrx,P_mtrx,t_inc,j); %solve for z
            SSM(k,j)=z(k,ODE_size); %only store the observable sensitivity result (GFP/MG) to SMM
        end

    end


    [SSM]=Normalize_SSM (DF,SSM,P,x,nstep);  % Identifiablity was estimated using an normalized SSM
    TF = isempty(IDP);
    if TF==0
        SSM(:,IDP(:,1))=0;
    end
end

