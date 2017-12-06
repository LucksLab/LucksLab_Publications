%Identifiability_N_Design_Main

%Construct 1: single AS no transcription
%Construct 2: single AS transcription
%Construct 3: double AS transcription
%Construct 4: mutt attn transcription
%Construct 5: H2-sTRSV-MG 
clear all
close all
CSTR_LV=1;
[ DF ] = Datafile(CSTR_LV );
R=DF.R;
IDP=[]; % Identifiable Parameters
fopen('Identifiable_Parm.txt','w');

%Loop to find identifiable param from CSTR_LV 1-5
while CSTR_LV <=R
    if CSTR_LV<=R && CSTR_LV>=1 
        DF = Datafile(CSTR_LV);%loads appropriate data presets for each construct level
    end
    [t,x] = Parameterization_Call_ODE(DF);%solves mass balance odes @ each construct level
    [SSM] = Parameterization_Numerical_SSM(DF,x,IDP); %this function uses the numerical method to find SSM
    [pset] = Identifiability (SSM); %array of of identifiable param @ each construct level
    IDP = [IDP;pset]; %append pset @ each cstr_lv to IDP

    print=['Construct Level ',num2str(CSTR_LV),'. Identified parameters: ',num2str(pset')]; 
    dlmwrite('Identifiable_Parm.txt',print,'delimiter','','-append'); %writes identifiable param @ each cstr_lv to file
    disp(print); 

    %Ploting
    A=abs(SSM);
    y=1:DF.Num_Parameters;
    x=0:5:100;
    figure(CSTR_LV) 
    imagesc(x,y,A')
    caxis([0 0.5])
    ylabel('Parameters'), xlabel('Time(min)')
    colorbar
    hold all
    filename=['CSTR_',num2str(CSTR_LV),'.fig'];
    savefig(filename);
    %Loop Counter
    CSTR_LV=CSTR_LV+1;
end


