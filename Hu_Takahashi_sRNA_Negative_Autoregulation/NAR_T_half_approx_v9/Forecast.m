function [tss_approx,norm_Data] = Forecast(norm_Data,storeSlope,k,slopeThresh,j,r)
    
    % This function approximates when the data set will reach steady state
    % if 5pts LSM fails.
    l=length(storeSlope);
    y = 100;        %initialize y before while loop
    tss_approx = k; %initialize steady state index before while loop
    
    % Create array of Xs
    index = [1:l];
  
    % Create Linear Fit for Second Derivative
    fitvars = polyfit(index(:), storeSlope(:),1);   %fit change in slope trend to linear line
    
    %polyfit returns array with [slope and intercept]
    m = fitvars(1); %slope
    b = fitvars(2); %intercept
    
    while y > slopeThresh
        y = m*tss_approx+b;        %approximate slope (y)
        if y < slopeThresh
            break;
        end
        tss_approx = tss_approx+1;  %add 1 as input back into while loop
    end
    
    % if projected steady state index out of data range, begin forecast
    if tss_approx > r
        for i = r+1:tss_approx
            norm_Data(i,j) = norm_Data(i-1,j)+(m*i+b);  %projects 1 point at a time to steady state index
        end
    end
end