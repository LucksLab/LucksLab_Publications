close all
clear all

CSTR_LV=1;
[DF] = Datafile_Analysis(CSTR_LV); %Load iteration number from struct
iter_num=DF.iter_num;

P = dlmread('P_solution.txt','\t');
P_size=size(P,1);

for i=1:P_size
    x=P(i,:);
    figure (i)
    hist(x,100);
    xlabel('value');
    ylabel('frequency');
    title_name=['Paremeter ',num2str(i)];
    title(title_name);
    hold all
end

