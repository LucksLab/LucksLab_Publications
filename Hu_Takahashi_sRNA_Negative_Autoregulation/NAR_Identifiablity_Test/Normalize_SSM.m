function [SSM]=Normalize_SSM (DF,SSM,P,x,nstep)

c=2;
NSSM=zeros(nstep./10,(DF.Num_Parameters));
ODE_size=DF.ODE_size;
f=DF.Timefactor;
for i=f+1:f:(nstep+1) % Pulling out coarse-grained time steps that correspond to experiments
    
    for j=1:(DF.Num_Parameters)
        NSSM(c,j)= SSM(i,j).*(P(j)./x(i,ODE_size)); %scaled by multiplying by pj/x to normalize
    end
    c=c+1;
end
SSM=NSSM;
end

