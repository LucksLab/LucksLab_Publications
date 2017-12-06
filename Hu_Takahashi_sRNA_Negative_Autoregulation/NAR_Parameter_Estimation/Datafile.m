%Datafile

function [ DF ] = Datafile(CSTR_LV,Identifiable_Param,P)

iter_num=1000;

t_i=0;
t_f= 7200; 
t_inc = 300; 
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

%initial guess of all parameters
P_i= dlmread('guess.txt','\t'); 
P_size=size(P_i,1);



%sets transcriptional rate of AS to be 8x of transcriptional rate of
%pre-cleaved RNA

if CSTR_LV==5  %CSTR_LV 5 has 2 mass balance odes
    ODE_size=2;
else %Other CSTR_LVs have 5 mass balance odes
    ODE_size=5;
end

%set initial condition to be 0 for all mass balance odes
IC=zeros(ODE_size,1);    
P_i_lb=P_i.*0.85;    %artificial range, within 15% of the best guess
P_i_ub=P_i.*1.15;   
a=P_i_lb;           %generate a set of P from 15% range about primary guess, first round, all random
b=P_i_ub ;
TF=isempty(P);
if TF==1 
    for j=1:(iter_num)
        for i=1:(P_size)
            P(j,i)=a(i)+(b(i)-a(i)).*rand(1,1);
        end
    end
end 

    
   
if CSTR_LV==1 %Single AS transcription rate = 0

    P(:,9)=0;
    P_i(9)=0;
    
else
    P(:,9)=8.*P(:,1);
    P_i(9)=8.*P_i(1);
end

 
DF.Initial_Parameters=P_i;

DF.Initial_Conditions=IC;
DF.Construct=CSTR_LV;
DF.Num_Parameters=P_size;
DF.t_i=t_i;
DF.t_f=t_f; 
DF.t_inc=t_inc; 
DF.nstep=nstep;               
DF.tspan=tspan;
DF.ODE_size=ODE_size;
DF.Parameter_library=P;
DF.iter_num=iter_num;
DF.Identifiable_Param=Identifiable_Param;

return;