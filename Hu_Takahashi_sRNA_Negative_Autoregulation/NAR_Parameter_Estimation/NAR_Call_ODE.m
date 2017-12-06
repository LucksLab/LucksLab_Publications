%Parameterization_Call_ODE
%This function uses ode15s to solve equations in MassbalanceEqns and
%returns the solution. 

function[t,x] = NAR_Call_ODE(DF)
P=DF.Initial_Parameters; 
tspan=DF.tspan;
CSTR_LV=DF.Construct;
ODE_size=DF.ODE_size;
x0=DF.Initial_Conditions;

% Solves Mass Balance Equations
[t,x]=ode15s(@(t,x) NAR_MassbalanceEqns(t,x,P,CSTR_LV,ODE_size),tspan, x0);
end

