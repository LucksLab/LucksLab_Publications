function [Exp_Data_DF] = Exp_Data_NAR( CSTR_LV,R )

    filename=['NAR_',num2str(CSTR_LV,R),'.txt'];
    Data = dlmread(filename, '\t'); 
    [r,c]=size(Data);
    time=Data(:,1);
    MG=Data(:,2:c);  
    stepsz=size(time,1);
    
    for i=1:(stepsz)
        avg(i)=sum(MG(i,:))/(c-1); 
        stdv(i)=std(MG(i,:)); 
    end
  
    Exp_Data_DF.Data=MG;
    Exp_Data_DF.avg=avg;
    Exp_Data_DF.stdv=stdv;
    Exp_Data_DF.timestep=stepsz;


end

