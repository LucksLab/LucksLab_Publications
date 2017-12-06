close all
clear all
hold on

%Enter filename
filename=strcat(input('\n Enter file name (no extensions): ','s'),'.txt');
data = dlmread(filename);

% promt user for time interval between each time point
t_int= input('\n Enter time interval (in min) between each data point:');

[r,c]=size(data);
scaling = max(max(data));       %Find Scaling factor

norm_Data = Norm_Data(data);    %Normalize data

%Storage Arrays
T_ss=[];         %Array to store steady state index
storeSlope=[];   %Array to store 5 points LSM slopes
    
%-----------------------CHANGE THIS-------------------------%
slopeThresh=0.005; %varies threshold value to call steady state
%-----------------------------------------------------------%

%------------------------SMOOTHING--------------------------%
wantSmooth = input('\n Need to filter noises in data points? (Y-1/N-0)');

    if wantSmooth == 1
        norm_Data(:,:) = Smoothing(norm_Data);
        wantSave = input('\n Want to save the filtered data in separate file? (Y-1/N-0)');
        if wantSave == 1
            filterName=strcat(input('\n Enter desired file name (no extensions): ','s'),'.txt');
            dlmwrite(filterName,norm_Data(:,:).*scaling,';'); %creates output file name with name desired
        end
    end
%-----------------------------------------------------------%

%sometimes it's nice to know whether your matlab is not responding or if it
%is doing actual work...
fprintf('\n begin analysis ... \n');

% Parse trials in columns (each column is a different data set)
for j=1:c

    storeSlope = zeros();   % refreshes 5pts LSM slopes for each trials

    for k=1:r-4     % Parse 5 time points in rows for each trials

        % Perform Linear Least Square Method on every 5 data points
        slope = LSM(norm_Data,k,j);
        storeSlope(k)=slope;    %store all slopes calculated for 2nd derivative test and forecasting

        % A slope within the set treshold is assumed to be in steady state
        if slope < slopeThresh
               T_ss(j)=(k);
               break;   %break loop move on to next column (data set)

        elseif k >= r-4 %apply second derivative test to approximate T_ss
            %get prepared for the most ugly function inputs ever...
            [forecast_tss,norm_Data] = Forecast(norm_Data,storeSlope,k,slopeThresh,j,r);
            T_ss(j) = forecast_tss;
        end
    end

    Half_SSval = norm_Data(T_ss(j),j)/2;    %calculate half steady state expression

    %loop through to find half steady state intervals
    for k = 1:r
        if norm_Data(k,j) > Half_SSval
            pt_high = k;
            pt_low = pt_high-1;
            break
        end
    end
    T_half(j) = double(Interpolate_Thalf(norm_Data(:,j),pt_low,pt_high,Half_SSval));   %interpolate between two points to get half time, set var type to double to get human number outputs...
end

    % time calc
    T_ss(:)=(T_ss(:)-1)*t_int;      % -1 because t = 0 min at index 1
    T_half(:)=(T_half(:)-1)*t_int;  % -1 because t = 0 min at index 1
    
    % output file
    outputname=strcat(input('\n Enter output file name (no extensions): ','s'),'.txt');
    dlmwrite(outputname,T_half(:),';'); %creates output file name with name desired

