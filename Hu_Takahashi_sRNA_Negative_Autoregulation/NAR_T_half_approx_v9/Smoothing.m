function [smooth_Data]= Smoothing(norm_Data)

%t= 0:5:120; %time step for tx-tl
%t= 0:20:320; %time step for invivo

[r,c]=size(norm_Data);
for i = 1:c 
    smooth_Data(:,i) = smooth(norm_Data(:,i));
    %plot(t,smooth_Data(:,i))  %time step above only relevant for plotting purposes
end
end