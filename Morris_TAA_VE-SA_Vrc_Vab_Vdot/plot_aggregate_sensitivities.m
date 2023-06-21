close all
sims=[0 1 2 3 4 5 6 7 8 9];
all_sens=[];
for n=1:length(sims)
    currentFile = sprintf('Sens_%d.mat',sims(n));
    load(currentFile)
    plot(sens_norm_col)
    all_sens=[all_sens sens_norm_col];
    hold on
end
medians=median(all_sens');
plot(medians,'k-')
set(gca,'XTick',1:numel(sens_norm_col),'XTickLabel',par_names)
xtickangle(60)
legend('0','1','2','3','4','5','6','7','8','9','median')

[sens_median,index] = sort(medians,'descend');
median_par_names=par_names(index,:)
f=figure;
for n=1:length(sims)
    sens_sort=all_sens(index,n)
    plot(sens_sort,'LineWidth',3);
    hold on
end
plot(sens_median,'k*-','LineWidth',3)
set(gca,'XTick',1:numel(medians),'XTickLabel',median_par_names)
xtickangle(60)
legend('Sim 0','Sim 1','Sim 2','Sim 3','Sim 4','Sim 5','Sim 6','Sim 7','Sim 8','Sim 9','median','Location','eastoutside')
%title('V_{rc},V_Aggregate Sensitivities','FontName','Arial','FontSize',30);
ylabel('Relative Sensitivity Index','FontSize',20);
set(gca,'FontSize',20);
f.Position=[200 200 1000 500];
saveas(gcf,'V_rc,V_ab,Vdot_Aggregate_Sensitivities.fig');
