close all
figure;
hold on
set(gca,'XTick',1:numel(Rsens),'XTickLabel',par_names);
ylabel('Relative Sensitivities');
Legend=cell(11,1);
iter=1;
for i=[0,1,3,5,6,7,9,10,11,12,13]
    simstr=num2str(i);
    str= ['Sens_' simstr '.mat'];
    load(str);
    newIsens=[];
    newRsens=[];
    for k=1:numel(Rsens)
    idx=find(Isens==k);
    newIsens(k)=Isens(idx);
    newRsens(k)=Rsens(idx);
    end
    plot(newRsens);
    set(gca,'FontSize',12);   
    Legend{iter}=strcat('Sim', simstr);
    iter=iter+1;
end
legend(Legend);