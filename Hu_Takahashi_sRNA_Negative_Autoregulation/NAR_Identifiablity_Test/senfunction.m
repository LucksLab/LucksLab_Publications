% the sensitivity coefficient matrix z(t), is a time-varying matrix that 
% encapsulates how sensitive xi is to a change in the pj

function [DZDT]=senfunction(t,z,J_mtrx,P_mtrx,t_inc,j) %M,PM)
 
        JM=J_mtrx(:,:); %gets Jacobian Matrix
        PM=P_mtrx(:,j); %gets GFP/MG row in P_Matrix
        
        DZDT=JM*z+PM; %forms ODE to solve for sensitivity coefficient z

end 