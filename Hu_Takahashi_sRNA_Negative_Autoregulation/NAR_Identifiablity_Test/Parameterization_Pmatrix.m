% df/dp

function [ P_mtrx ] = Parameterization_Pmatrix(P,x,k,CSTR_LV,ODE_size,P_size,j)

P_mtrx=zeros(ODE_size,P_size);     % initialize the P_matrix
t=0;
P_holder=P;
x=x(k,:); %get x at time point k
for i=1:ODE_size
        P=P_holder;
        F=zeros(4,1);
        h=P(j)*0.01;
        %Gets O(4) central difference on dfi/dpj
        if h~=0
            P(j)=P_holder(j)+2*h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(1)=f(i);
            P(j)=P_holder(j)+h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(2)=f(i); 
            P(j)=P_holder(j)-h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(3)=f(i);
            P(j)=P_holder(j)-2*h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(4)=f(i);
            %Store approx. dfi/dpj into P_matrix
            P_mtrx(i,j)= (-F(1)+8.*F(2)-8.*F(3)+F(4))./(12.*h);   
        end

end


end

