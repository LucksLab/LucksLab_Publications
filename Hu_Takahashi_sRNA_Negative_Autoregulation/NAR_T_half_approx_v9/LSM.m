% Least Square Method

function [LSM_Slope] = LSM(norm_data,k,j)

    %least square on 5 data points
    n=5;
    t=k:k+4;
    GFP = [norm_data(k,j),norm_data(k+1,j),norm_data(k+2,j),norm_data(k+3,j),norm_data(k+4,j)];
    
    St=sum(t);
    SGFP=sum(GFP);
    StGFP=sum(t.*GFP);
    Stt=sum(t.*t);
    
    %calculate slope using least square formula
    slope=abs((n.*StGFP-St.*SGFP)/(n.*Stt-St.*St));
    %b = (SGFP/n) - (slope*St/n);
    
    LSM_Slope = slope;
end