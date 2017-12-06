function [Exp_Data_DF] = Exp_Data( CSTR_LV,R )

    for i=1:R
        filename=['CSTR_LV_',num2str(i),'.txt'];
        Data = dlmread(filename, '\t'); 
        GFP_max(i)=max(max(Data));
        GFP_min(i)=min(min(Data));
    end

    filename=['CSTR_LV_',num2str(CSTR_LV,R),'.txt'];
    Data = dlmread(filename, '\t'); 
    [r,c]=size(Data);
    time=Data(:,1);
    GFP=Data(:,2:c);  
    stepsz=size(time,1);
    
    for i=1:(stepsz)
        Upbound(i)=max(GFP(i,:));
        lowbound(i)=min(GFP(i,:));
        avg(i)=sum(GFP(i,:))/(c-1); 
    end

    Exp_Data_DF.Data=GFP;
    Exp_Data_DF.avg=avg';
    Exp_Data_DF.Upbound=Upbound';
    Exp_Data_DF.Lowbound=lowbound';
    Exp_Data_DF.timestep=stepsz;
    Exp_Data_DF.Max=max(GFP_max);
    Exp_Data_DF.Min=min(GFP_min);

end

