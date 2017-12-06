%Interpolate T1/2

function [t_half] = Interpolate_Thalf(array,pt_low,pt_high,Half_SSval)
    
    % interpolate 2 pts
    x = pt_low:pt_high;
    y = array(pt_low:pt_high);
    
    fitvars = polyfit(x(:), y(:),1);    %fit linear line to 2 point
    m = fitvars(1);
    b = fitvars(2);
    
    syms x
    t_half = solve(m*x+b == Half_SSval, x); %solve for half steady state value
end