%Datafile

function [ DF ] = Datafile(CSTR_LV )

t_i=0;
t_f= 6000; 
t_inc = 30; 
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

exp_t_inc=300; % Time interval of the trajectory in the experiments, in seccond
f=exp_t_inc./t_inc;

%initial guess of all parameters
P_i= dlmread('guess.txt','\t'); 
P_size=size(P_i,1);
 
R=5; %R is your total number of constructs

%sets transcriptional rate of AS to be 8x of transcriptional rate of
%pre-cleaved RNA
P_i(9)=8*P_i(1);

if CSTR_LV==1 %Single AS transcription rate = 0
    P_i(9)=0;
end

if CSTR_LV==5  %CSTR_LV 5 has 2 mass balance odes
    ODE_size=2;
else %Other CSTR_LVs have 5 mass balance odes
    ODE_size=5;
end

%set initial condition to be 0 for all mass balance odes
IC=zeros(ODE_size,1);    

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
DF.R=R;
DF.Timefactor=f;
return;