%df/dx

function [ J_mtrx ] = Parameterization_Jacobian(P,x,k,CSTR_LV,ODE_size)

J_mtrx=zeros(ODE_size,ODE_size);     % initialize the Jacobian matrix
t=0;
X=x; % Use X as x storage

for i=1:ODE_size
    for j=1:ODE_size
        F=zeros(4,1);
        h=X(k,j).*0.01;       
        if h~=0
            %Gets x at time point k
            x=X(k,:);
            %Gets O(4) central difference on dfj/dxi
            x(j)=X(k,j)+2*h;
        	[f]= MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(1)=f(i);
            x(j)=X(k,j)+h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(2)=f(i); 
            x(j)=X(k,j)-h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(3)=f(i);
            x(j)=X(k,j)-2*h;
            [f] = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size);
            F(4)=f(i);
            %Stores appox. dfj/dxi into J matrix
            J_mtrx(i,j)= (-F(1)+8*F(2)-8*F(3)+F(4))./(12*h);   
        end
    end
end
end

