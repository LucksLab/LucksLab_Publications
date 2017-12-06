    %This function contains the massbalance ODEs of five parameterization
    %experiments
function f = MassbalanceEqns(t,x,P,CSTR_LV,ODE_size)

f = zeros(ODE_size,1);

        if CSTR_LV==1 || CSTR_LV==2
            %single antisense
            f(1)=P(9)-P(5).*x(1)-P(6).*x(1);
            f(2)=P(6).*x(1)-P(5).*x(2);
            f(3)=P(1).*(1-P(16)).*(1-(x(2)./(P(2)+x(2))))-P(12).*x(3);
            f(4)=P(10).*x(3)-P(11).*x(4);
            f(5)=P(11).*x(4);
            
        elseif CSTR_LV==3
            %double antisense
           f(1)=P(9)-P(15).*x(1)-P(7).*x(1);
            f(2)=P(7).*x(1)-P(15).*x(2);
            f(3)=P(1).*(1-P(16)).*(1-(x(2)./(P(3)+x(2))))-P(12).*x(3);
            f(4)=P(10).*x(3)-P(11).*x(4);
            f(5)=P(11).*x(4);     
            
        elseif CSTR_LV==4
            %mutt attn antisense
            f(1)=P(9)-P(5).*x(1)-P(6).*x(1);
            f(2)=P(6).*x(1)-P(5).*x(2);
            f(3)=P(1).*(1-P(16)).*(1-(x(2)./(P(4)+x(2))))-P(12).*x(3);
            f(4)=P(10).*x(3)-P(11).*x(4);
            f(5)=P(11).*x(4);    
            
        elseif CSTR_LV==5
            %MG
            f(1)=P(13).*(1-P(16))-P(14).*x(1)-P(8).*x(1);
            f(2)=P(8).*x(1)-P(14).*x(2);
            
        end
        
        
 end
