function [DF] = Datafile_Prediction(CSTR_LV)

iter_num=100;
Total_constructs=4;
t_i=0;
t_f= 7200; 
t_inc = 300; 
nstep = (t_f-t_i)/t_inc;               
tspan = t_i:t_inc:t_f;

P = dlmread('P_solution.txt','\t');
P_size=size(P,1);

ODE_size=4;
IC=zeros(ODE_size,1);  

DF.Parameter_library=P;
DF.Initial_Conditions=IC;
DF.Construct=CSTR_LV;
DF.Num_Parameters=P_size;
DF.t_i=t_i;
DF.t_f=t_f; 
DF.t_inc=t_inc; 
DF.nstep=nstep;               
DF.tspan=tspan;
DF.ODE_size=ODE_size;

DF.iter_num=iter_num;
DF.Total_constructs=Total_constructs;
end


