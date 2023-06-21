addpath('C:/Users/rfita/OneDrive/Documents/MATLAB/Chest Vest Code/Copy_of_Sensitivity_TAA_AJP_Fixed_Sampling/Morris_VE-SA_Vdot/');
fig=figure(1);

pl=tiledlayout(1,3);
pl.Padding="compact";
pl.TileSpacing="tight";

colormatrix=zeros(10,3);
colormatrix(1,:)=[0.5 0.5 0.5];
colormatrix(2,:)=[0.3 0.7 0.7];
colormatrix(3,:)=[0.1 0.9 0.9];
colormatrix(4,:)=[0.3 1 0.5];
colormatrix(5,:)=[0.5 0.7 0.3];
colormatrix(6,:)=[0.7 0.5 0.1];
colormatrix(7,:)=[0.9 0.3 0.3];
colormatrix(8,:)=[1 0.1 0.5];
colormatrix(9,:)=[0.9 0.3 0.7];
colormatrix(10,:)=[0.7 0.5 0.9];

nexttile
sims=[0 1 2 3 4 5 6 7 8 9];
all_sens=[];
for n=1:length(sims)
    currentFile = sprintf('Sens_%d.mat',sims(n));
    load(currentFile)
    all_sens=[all_sens sens_norm_col];
end
medians=median(all_sens');
[sens_median,index] = sort(medians,'descend');
median_par_names=par_names(index,:)
for n=1:length(sims)
    sens_sort=all_sens(index,n)
    plot(sens_sort,'LineWidth',1,"Color",colormatrix(n,:));
    hold on
end
plot(sens_median,'k*-','LineWidth',3)
hold off
set(gca,'XTick',1:numel(medians),'XTickLabel',median_par_names)
axis([1 7 0 15]);
xtickangle(-30)
ylabel('Sensitivity Index','FontSize',20);
set(gca,'FontSize',20);
title('$\dot{V}$','Interpreter','latex','FontWeight','bold');
tex=annotation('textbox','String','A');
tex.Parent=fig.CurrentAxes;
tex.Position=[0.5 22 3 2];
tex.FontSize=20;
tex.FontWeight='bold';
tex.LineStyle="none";


rmpath('C:/Users/rfita/OneDrive/Documents/MATLAB/Chest Vest Code/Copy_Of_Sensitivity_TAA_AJP_Fixed_Sampling/Morris_VE-SA_Vdot/');
addpath('C:\Users\rfita\OneDrive\Documents\MATLAB\Chest Vest Code\Copy_Of_Sensitivity_TAA_AJP_Fixed_Sampling\Morris_TAA-SA_Vrc_Vab\');
nexttile
sims=[0 1 2 3 4 5 6 7 8 9];
all_sens=[];
for n=1:length(sims)
    currentFile = sprintf('Sens_%d.mat',sims(n));
    load(currentFile)
    all_sens=[all_sens sens_norm_col];
end
medians=median(all_sens');
[sens_median,index] = sort(medians,'descend');
median_par_names=par_names(index,:)
for n=1:length(sims)
    sens_sort=all_sens(index,n)
    plot(sens_sort,'LineWidth',1,'Color',colormatrix(n,:));
    hold on
end
plot(sens_median,'k*-','LineWidth',3)
set(gca,'XTick',1:numel(medians),'XTickLabel',median_par_names)
xtickangle(-30)
axis([1 7 0 15]);
set(gca,'FontSize',20);
title('$V_{rc}$ and $V_{ab}$','Interpreter','latex','FontWeight','bold');
tex=annotation('textbox','String','B');
tex.Parent=fig.CurrentAxes;
tex.Position=[0.5 22 3 2];
tex.FontSize=20;
tex.FontWeight='bold';
tex.LineStyle="none";



rmpath('C:\Users\rfita\OneDrive\Documents\MATLAB\Chest Vest Code\Copy_Of_Sensitivity_TAA_AJP_Fixed_Sampling\Morris_TAA-SA_Vrc_Vab\');
addpath('C:\Users\rfita\OneDrive\Documents\MATLAB\Chest Vest Code\Copy_Of_Sensitivity_TAA_AJP_Fixed_Sampling\Morris_TAA_VE-SA_Vrc_Vab_Vdot\');
nexttile
sims=[0 1 2 3 4 5 6 7 8 9];
all_sens=[];
for n=1:length(sims)
    currentFile = sprintf('Sens_%d.mat',sims(n));
    load(currentFile)
    all_sens=[all_sens sens_norm_col];
end
medians=median(all_sens');

[sens_median,index] = sort(medians,'descend');
median_par_names=par_names(index,:)
for n=1:length(sims)
    sens_sort=all_sens(index,n)
    plot(sens_sort,'LineWidth',1,'Color',colormatrix(n,:));
    hold on
end
plot(sens_median,'k*-','LineWidth',3)
set(gca,'XTick',1:numel(medians),'XTickLabel',median_par_names)
axis([1 8 0 15]);
legend('Sim 0','Sim 1','Sim 2','Sim 3','Sim 4','Sim 5','Sim 6','Sim 7','Sim 8','Sim 9','median','Location','eastoutside')
xtickangle(-30)
set(gca,'FontSize',20);
title('$\dot{V},V_{rc}$, and $V_{ab}$','Interpreter','latex');
tex=annotation('textbox','String','C');
tex.Parent=fig.CurrentAxes;
tex.Position=[0.5 22 3 2];
tex.FontSize=20;
tex.FontWeight='bold';
tex.LineStyle="none";
fig.Position=[200 200 1200 400];