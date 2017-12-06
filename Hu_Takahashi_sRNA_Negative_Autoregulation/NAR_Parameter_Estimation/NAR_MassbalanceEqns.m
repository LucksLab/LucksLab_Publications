    %This function contains the massbalance ODEs of five parameterization
    %experiments
function f = NAR_MassbalanceEqns(t,x,P,CSTR_LV,ODE_size)

f = zeros(ODE_size,1);

        if CSTR_LV==1 
            %single antisense NAR
            f(1)=P(13).*(1-P(16)).*(1-(x(2)./(P(2)+x(2))))-P(5).*x(1)-P(6).*x(1);
            f(2)=P(6).*x(1)-P(5).*x(2);
            f(3)=P(13).*(1-P(16)).^2.*(1-(x(2)./(P(2)+x(2))))-P(14).*x(3)-P(8).*x(3);
            f(4)=P(8).*x(3)-P(14).*x(4);
        elseif CSTR_LV==2
            %single antisense NAR CTRL
            f(1)=P(13).*(1-P(16)).*(1-(x(2)./(P(4)+x(2))))-P(5).*x(1)-P(6).*x(1);
            f(2)=P(6).*x(1)-P(5).*x(2);
            f(3)=P(13).*(1-P(16)).^2.*(1-(x(2)./(P(4)+x(2))))-P(14).*x(3)-P(8).*x(3);
            f(4)=P(8).*x(3)-P(14).*x(4);
            
        elseif CSTR_LV==3
            %double antisense NAR
            f(1)=P(13).*(1-P(16)).*(1-(x(2)./(P(3)+x(2))))-P(15).*x(1)-P(7).*x(1);
            f(2)=P(7).*x(1)-P(15).*x(2);
            f(3)=P(13).*(1-P(16)).^3.*(1-(x(2)./(P(3)+x(2))))-P(14).*x(3)-P(8).*x(3);
            f(4)=P(8).*x(3)-P(14).*x(4);
            
        elseif CSTR_LV==4
             %double antisense NAR CTRL
            f(1)=P(13).*(1-P(16)).*(1-(x(2)./((P(3)./P(2).*P(4))+x(2))))-P(15).*x(1)-P(7).*x(1);
            f(2)=P(7).*x(1)-P(15).*x(2);
            f(3)=P(13).*(1-P(16)).^3.*(1-(x(2)./((P(3)./P(2).*P(4))+x(2))))-P(14).*x(3)-P(8).*x(3);
            f(4)=P(8).*x(3)-P(14).*x(4);
            
        end
        
 end
