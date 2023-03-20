function P_FRC=test_curves_preterm_freeRCAB(curvepars)

%Call Parameters
TLC=curvepars(1);
RV=curvepars(2);
plopn=curvepars(3);
plran=curvepars(4);
beta=curvepars(5);
gamma=curvepars(6);
k=curvepars(7);
vrcstr=curvepars(8);
vabstr=curvepars(9);
prcMult=curvepars(10);

VC=(TLC-RV);

Pel=transpose(-25:.01:35);
Pcw=transpose(-25:.01:35);

%Rib Cage Auxiliary Parameters
prcsl=22.6974/log((exp(0.008/vrcstr)-1)/(exp(0.0001/vrcstr)-1))*prcMult; 
prctra=-prcsl*log(exp(0.008/vrcstr)-1);

%Rib Cage Volume
frc=vrcstr*log(1+exp((Pcw-prctra)/prcsl))/VC;
Vrc=frc*VC;

%Abdomen Auxiliary Parameters
pabsl=5.6/log((exp(0.00476496/vabstr)-1)/(exp(0.00116/vabstr)-1)); %2.92
pabtra=-pabsl*log(exp(0.00116/vabstr)-1); %3.18

%Abdomen Volume
fab=vabstr*log(1+exp((Pcw-pabtra)/pabsl))/VC;
Vab=fab*VC;

%Alveoli/Lung
alpha=((1+exp(plopn/plran))*beta-gamma)/exp(plopn/plran);
Frec=alpha+(gamma-alpha)./(1+exp(-(Pel-plopn)/plran));
Vel=VC*(1-exp(-k*Pel));
VA=Frec.*Vel+RV;
Vcw=Vrc+Vab+RV;

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

index=find(abs(Ptot)<0.03,1,'last');
FRC=VolN(index);
P_FRC=PelN(index);

Pel_range=Pel;
Vcw_range=Vcw;
VA_range=VA;
save PVcurves.mat Pel_range Vcw_range VA_range Ptot VolN FRC P_FRC