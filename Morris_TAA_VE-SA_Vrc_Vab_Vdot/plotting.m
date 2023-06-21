% clear all
close all

set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesColorOrder',[  0 0 0; 0    0.4470    0.7410;0.8500    0.3250    0.0980;    ])
set(0,'defaultAxesLineStyleOrder',{'-',':'}) 

%         0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840
    
% %%
statesCell=cell(6,6);

load states_highCw_normRu_zeroPao.mat
statesCell(1,1)={Vdot*1000};
statesCell(1,2)={VA*1000};
statesCell(1,3)={PA};
statesCell(1,4)={Pldyn};
statesCell(1,5)={Ppl};
statesCell(1,6)={t};

load states_highCw_highRu_zeroPao.mat
statesCell(2,1)={Vdot*1000};
statesCell(2,2)={VA*1000};
statesCell(2,3)={PA};
statesCell(2,4)={Pldyn};
statesCell(2,5)={Ppl};% suff='highCw';
statesCell(2,6)={t};

load states_highCw_normRu_cpapPao.mat
statesCell(3,1)={Vdot*1000};
statesCell(3,2)={VA*1000};
statesCell(3,3)={PA};
statesCell(3,4)={Pldyn};
statesCell(3,5)={Ppl};% suff='highCw';
statesCell(3,6)={t};

load states_lowCw_normRu_zeroPao.mat
statesCell(4,1)={Vdot*1000};
statesCell(4,2)={VA*1000};
statesCell(4,3)={PA};
statesCell(4,4)={Pldyn};
statesCell(4,5)={Ppl};% suff='highCw';
statesCell(4,6)={t};
%
load states_lowCw_highRu_zeroPao.mat
statesCell(5,1)={Vdot*1000};
statesCell(5,2)={VA*1000};
statesCell(5,3)={PA};
statesCell(5,4)={Pldyn};
statesCell(5,5)={Ppl};% suff='highCw';
statesCell(5,6)={t};

load states_lowCw_normRu_cpapPao.mat
statesCell(6,1)={Vdot*1000};
statesCell(6,2)={VA*1000};
statesCell(6,3)={PA};
statesCell(6,4)={Pldyn};
statesCell(6,5)={Ppl};% suff='highCw';
statesCell(6,6)={t};


time=statesCell{1,6}(101:601);

figure
plot(time,statesCell{1,1}(101:601),time,statesCell{2,1}(101:601),time,statesCell{3,1}(101:601),time,statesCell{4,1}(101:601),time,statesCell{5,1}(101:601),time,statesCell{6,1}(101:601))
title('Airflow')
% xlabel('Time, s')
ylabel('dV/dt, ml/s')
legend('High C_w','High C_w + high R_u','High C_w + CPAP','Low C_w','Low C_w + high R_u','Low C_w + CPAP','Location','best')
saveas(gcf, 'Vdot', 'fig')

figure
plot(time,statesCell{1,2}(101:601),time,statesCell{2,2}(101:601),time,statesCell{3,2}(101:601),time,statesCell{4,2}(101:601),time,statesCell{5,2}(101:601),time,statesCell{6,2}(101:601))
title('Alveolar volume')
% xlabel('Time, s')
ylabel('V_A, ml')
legend('High C_w','High C_w + high R_u','High C_w + CPAP','Low C_w','Low C_w + high R_u','Low C_w + CPAP','Location','best')
saveas(gcf, 'VA', 'fig')

figure
plot(time,statesCell{1,5}(101:601),time,statesCell{2,5}(101:601),time,statesCell{3,5}(101:601),time,statesCell{4,5}(101:601),time,statesCell{5,5}(101:601),time,statesCell{6,5}(101:601))
title('Pleural pressure')
% xlabel('Time, s')
ylabel('P_{pl}, cm H_2O')
legend('High C_w','High C_w + high R_u','High C_w + CPAP','Low C_w','Low C_w + high R_u','Low C_w + CPAP','Location','best')
saveas(gcf, 'Ppl', 'fig')

figure
plot(time,statesCell{1,3}(101:601),time,statesCell{2,3}(101:601),time,statesCell{3,3}(101:601),time,statesCell{4,3}(101:601),time,statesCell{5,3}(101:601),time,statesCell{6,3}(101:601))
title('Alveolar pressure')
% xlabel('Time, s')
ylabel('P_A, cm H_2O')
legend('High C_w','High C_w + high R_u','High C_w + CPAP','Low C_w','Low C_w + high R_u','Low C_w + CPAP','Location','best')
saveas(gcf, 'PA', 'fig')

figure
plot(time,statesCell{1,4}(101:601),time,statesCell{2,4}(101:601),time,statesCell{3,4}(101:601),time,statesCell{4,4}(101:601),time,statesCell{5,4}(101:601),time,statesCell{6,4}(101:601))
title('Dynamic lung elastic recoil')
% xlabel('Time, s')
ylabel('P_{l,dyn}, cm H_2O')
legend('High C_w','High C_w + high R_u','High C_w + CPAP','Low C_w','Low C_w + high R_u','Low C_w + CPAP','Location','best')
saveas(gcf, 'Pel', 'fig')
return

%%
load highCw.mat
Cdynhigh=Cdynsave; VThigh=VTsave;

% figure
% set(gcf,'defaultAxesColorOrder',[ 0.8500    0.3250    0.0980; .5 0 0])
% subplot(1,2,1)
% yyaxis left
% plot([10:length(VAminsave)]'/3600,VAminsave(10:end)*1000)
% title('Lung volumes, high C_w no intervention (S1)')
% xlabel('Time, h')
% ylabel('End-expiratory lung volume, ml')
% axis([-Inf Inf 23 30])
% yyaxis right
% plot([10:length(VTsave)]'/3600,VTsave(10:end)*1000); %#ok<NBRAK>
% ylabel('Tidal volume V_T, ml')
% legend('EELV','V_T')
% 
load lowCw.mat
Cdynlow=Cdynsave; VTlow=VTsave;

% set(gcf,'defaultAxesColorOrder',[  0    0.4470    0.7410; 0 0 .5])
% subplot(1,2,2)
% yyaxis left
% plot([10:length(VAminsave)]'/3600,VAminsave(10:end)*1000)
% title('Lung volumes, low C_w no intervention (S7)')
% xlabel('Time, h')
% ylabel('End-expiratory lung volume, ml')
% axis([-Inf Inf 23 30])
% yyaxis right
% plot([10:length(VTsave)]'/3600,VTsave(10:end)*1000); %#ok<NBRAK>
% ylabel('Tidal volume V_T, ml')
% legend('EELV','V_T')
% 
% saveas(gcf, 'lungvolumesCw', 'fig')

load highCw_Ru.mat;CdynRu=Cdynsave; VTRu=VTsave;
load highCw_CPAP.mat;CdynCPAP=Cdynsave; VTCPAP=VTsave;
load highCw_CPAP5.mat;CdynCPAP5=Cdynsave; VTCPAP5=VTsave;
load highCw_CPAP3.mat;CdynCPAP3=Cdynsave; VTCPAP3=VTsave;

% figure
% plot([10:length(Cdynhigh)]'/3600,Cdynhigh(10:end)*1000,[10:length(Cdynlow)]'/3600,Cdynlow(10:end)*1000,[10:length(CdynCPAP)]'/3600,CdynCPAP(10:end)*1000,...
%     [10:length(CdynCPAP5)]'/3600,CdynCPAP5(10:end)*1000,[10:length(CdynCPAP3)]'/3600,CdynCPAP3(10:end)*1000); %#ok<NBRAK>
% title('Dynamic lung compliance')
% xlabel('Time, h')
% ylabel('C_L, ml / cm H_2O')
% legend('High Cw','Low Cw','High Cw + CPAP', 'High Cw + CPAP5', 'High Cw + CPAP3')

figure
plot([10:length(VThigh)]'/3600,VThigh(10:end)*1000,[10:length(VTlow)]'/3600,VTlow(10:end)*1000,[10:length(VTCPAP)]'/3600,VTCPAP(10:end)*1000,...
  [10:length(VTCPAP5)]'/3600,VTCPAP5(10:end)*1000,[10:length(VTCPAP3)]'/3600,VTCPAP3(10:end)*1000); %#ok<NBRAK>
title('Tidal volume')
xlabel('Time, h')
ylabel(' V_T, ml')
 legend('High Cw','Low Cw','High Cw + CPAP', 'High Cw + CPAP5', 'High Cw + CPAP3')
