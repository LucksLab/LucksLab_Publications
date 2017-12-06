% Normalize Experimental Data

function [Norm_Data] = Norm_Data(data)

    GFP_Max = max(max(data));   %find max
    Norm_Data = data./GFP_Max;  %normalize
    
end