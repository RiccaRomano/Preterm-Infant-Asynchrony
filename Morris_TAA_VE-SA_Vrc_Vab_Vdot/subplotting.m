close all

%% Figure 2

% pvhf=hgload('PVcurvesHighCw.fig');
% pvlf=hgload('PVcurvesLowCw.fig');
% 
% figure
% h(1)=subplot(1,2,1);
% axis([-10 35 RV*1000 TLC0*1000])
% h(2)=subplot(1,2,2);
% axis([-10 35 RV*1000 TLC0*1000])
% 
% copyobj(allchild(get(pvhf,'CurrentAxes')),h(1));
% copyobj(allchild(get(pvlf,'CurrentAxes')),h(2));
% 
% t(1)=title(h(1),'High C_w');
% t(2)=title(h(2),'Low C_w');
% 
% y(1)=ylabel(h(1),'Volume, ml');
% % y(2)=ylabel(h(2),'Volume, ml');
% 
% x(2)=xlabel(h(2),'Pressure, cm H_2O');
% l(2)=legend(h(2),'Chest wall','Lung','Total respiratory','Tidal loop','FRC','Location','best');
%% Figure 3
VAf=hgload('VA.fig');
Vdotf=hgload('Vdot.fig');
PAf=hgload('PA.fig');
Pelf=hgload('Pel.fig');
Pplf=hgload('Ppl.fig');

figure
h(1)=subplot(3,2,1);
h(2)=subplot(3,2,2);
h(3) = subplot(3,2,3);
h(4) = subplot(3,2,4);
h(5) = subplot(3,2,5); % the last (odd) axes

copyobj(allchild(get(Vdotf,'CurrentAxes')),h(1));
copyobj(allchild(get(Pplf,'CurrentAxes')),h(2));
copyobj(allchild(get(VAf,'CurrentAxes')),h(3));
copyobj(allchild(get(Pelf,'CurrentAxes')),h(4));
copyobj(allchild(get(PAf,'CurrentAxes')),h(5));
return
t(1)=title(h(1),'Airflow');
t(2)=title(h(2),'Pleural pressure');
t(3)=title(h(3),'Alveolar volume');
t(4)=title(h(4),'Dynamic elastic lung recoil');
t(5)=title(h(5),'Alveolar pressure');

y(1)=ylabel(h(1),'dV/dt, ml/s');
y(2)=ylabel(h(2),'P_{pl}, cm H_2O');
y(3)=ylabel(h(3),'V_A, ml');
y(4)=ylabel(h(4),'P_{l,dyn}, cm H_2O');
y(5)=ylabel(h(5),'P_A, cm H_2O');

x(5)=xlabel(h(5),'Time, s');
% l(5)=legend(h(5),'High C_w','Low C_w','High C_w + high R_u','Low C_w + high R_u','High C_w + CPAP','Low C_w + CPAP','Location','best')
l(5)=legend(h(5),'High C_w','High C_w + high R_u','High C_w + CPAP','Low C_w','Low C_w + high R_u','Low C_w + CPAP','Location','best');


pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)])

%% Figure 5

cdynf=hgload('Cdyn.fig');
vtf=hgload('VTdyn.fig');

figure
h(1)=subplot(1,2,1);
h(2)=subplot(1,2,2);

copyobj(allchild(get(cdynf,'CurrentAxes')),h(1));
copyobj(allchild(get(vtf,'CurrentAxes')),h(2));

t(1)=title(h(1),'Dynamic lung compliance');
t(2)=title(h(2),'Tidal volume');

y(1)=ylabel(h(1),'C_L, ml / cm H_2O');
y(2)=ylabel(h(2),'V_T, ml');

x(2)=xlabel(h(2),'Time, h');
l(2)=legend(h(2),'High Cw','Low Cw','High Cw + CPAP', 'High Cw + CPAP5', 'High Cw + CPAP3','Location','best');