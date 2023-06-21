function P_FRC=test_curves_preterm_freeRCAB(curvepars)
%close all
global drcMult RV TLC VC Phalf Ptau k
global beta gamma alpha
global brc crc drc bab cab dab

TLC=curvepars(1);
RV=curvepars(2);
Phalf=curvepars(3);
Ptau=curvepars(4);
beta=curvepars(5);
gamma=curvepars(6);
%alpha=curvepars(7);
k=curvepars(7);
brc=curvepars(8);
bab=curvepars(9);
drcMult=curvepars(10);

VC=(TLC-RV);
Pel=[-20:.01:35]';
Pcw=[-20:.01:35]';

%%% RC
% brc=0.00425; %normal bw=0.00425
drc=22.6974/log((exp(0.008/brc)-1)/(exp(0.0001/brc)-1))*drcMult; 
crc=-drc*log(exp(0.008/brc)-1);%

frc=brc*log(1+exp((Pcw-crc)/drc))/VC;
Vrc=frc*VC;
Vrcrv=Vrc+RV;

%%% AB
% bab=0.004;
dab=5.6/log((exp(0.00476496/bab)-1)/(exp(0.00116/bab)-1))*1; %2.92
cab=-dab*log(exp(0.00116/bab)-1); %3.18

fab=bab*log(1+exp((Pcw-cab)/dab))/VC;
Vab=fab*VC;
Vabrv=Vab+RV;

% figure
% plot(Pcw,Vabrv,Pcw,Vrcrv)
% axis([-20 20 RV TLC])
% line([0 0],[0 TLC],'Color',[.7 .7 .7])

alpha=((1+exp(Phalf/Ptau))*beta-gamma)/exp(Phalf/Ptau);
Frec=alpha+(gamma-alpha)./(1+exp(-(Pel-Phalf)/Ptau));
Vel=VC*(1-exp(-k*Pel));
VA=Frec.*Vel+RV;
Vcw=Vrc+Vab+RV; %aw+bw*log(1+exp((Pcw-cw)/dw));

itemp_start=find(Pcw==-15);
iVcw_start=find(abs(Vcw-Vcw(itemp_start))<0.0005, 1, 'last' ); %volume Vcw at that index
iVA_start=find(abs(VA-Vcw(itemp_start))<.0005, 1, 'last' ); %volume VA at the same index

itemp_end=find(abs(Vcw-TLC)<0.02,1, 'last');
iVcw_end=find(abs(Vcw-Vcw(itemp_end))<0.01, 1, 'last' );
iVA_end=find(abs(VA-Vcw(itemp_end))<1, 1, 'last' );

pcwt=Pcw(iVcw_start:iVcw_end);
vcwt=Vcw(iVcw_start:iVcw_end);
pelt=Pel(iVA_start:iVA_end);
vat=VA(iVA_start:iVA_end);

VolN=linspace(min(vat),max(vat),5000);
PelN=interp1(vat,pelt,VolN);
PcwN=interp1(vcwt,pcwt,VolN);
Ptot=PelN+PcwN;

index=find(abs(Ptot)<0.03,1,'last'); %0.01 or 0.02
FRC=VolN(index);
P_FRC=PelN(index);

Pel_range=Pel;Vcw_range=Vcw;VA_range=VA;
save PVcurves.mat Pel_range Vcw_range VA_range Ptot VolN FRC P_FRC
% return
% figure
% plot(Pcw,Vcw,'b',Pel(2001:end),VA(2001:end),'r',Ptot,VolN,'m');%Pel,Vel+RV,'g',
% hold on
% plot(PcwN(index),FRC,'k*',PelN(index),FRC,'k*')
% axis([-10 35 RV TLC])
% line([0 0],[0 TLC],'Color',[.7 .7 .7])
% line([-10 35],[RV RV],'Color',[.4 .4 .4],'Linestyle','--')
% line([-10 35],[FRC FRC],'Color',[.4 .4 .4],'Linestyle',':')
% xlabel('Pressure')
% ylabel('Volume, ml')
% legend('Chest wall','Lung','Total (chest plus lung)','FRC','Location','best')


%TLC=.001*63; 
%RV=.001*23; 
% alpha_cw=0.25;%
% volume0Pres=(alpha_cw*VC+RV);
