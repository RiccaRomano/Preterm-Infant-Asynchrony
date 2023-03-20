function [dpdt] =Model(t,p,T,tprev,pars)

%Call Parameters
TLC=pars(1);
RV=pars(2);
plopn=pars(3);
plran=pars(4);
beta=pars(5);
gamma=pars(6);
k=pars(7);
vrcstr=pars(8);
vabstr=pars(9);
prcMult=pars(10);
f=pars(12);
T=pars(13);
Pfrac=pars(14);
Apic=pars(15);
Amus=pars(16);
Rsm=pars(17);
Rsd=pars(18);
Ks=pars(19);
Rum=pars(20);
Ku=pars(21);
I=pars(22); 
pcmax=pars(23);
pcran=pars(24);
Kc=pars(25);
Vcmax=pars(26);
Cve=pars(27);
Rve=pars(28);
Rrc=pars(29);
Rab=pars(30);
Pao=pars(31);

%Call ODEs
Vdot = p(1);
Pel = p(2);
Vc = p(3);
Pve = p(4);
Vrc = p(5);

Pldyn=Pel+Pve; %Dynamic Lung Compliance
Ptm = pcmax-pcran*log(Vcmax./(Vc)-1); %Transmural Pressure

% Driving Pressures
Pmus = (Amus*cos(2*pi*f*(t-tprev)) + -Amus); %Diaphragm Pressure
Pdi_co=Pfrac*Pmus; %Costal Diaphragm
Pdi_cr=Pmus-Pdi_co; %Crural Diaphragm
Pic = (Apic*cos(2*pi*f*(t-tprev)) + -Apic); %Intercostal Pressure

%Lung Compartment
VC  = TLC-RV;
alpha=((1+exp(plopn/plran))*beta-gamma)/exp(plopn/plran);
Frec=alpha+(gamma-alpha)./(1+exp(-(Pel-plopn)/plran));
Vel=VC*(1-exp(-k*Pel));
%Alveolar Volume
VA=Frec*Vel+RV;

%Chest Wall and Abdomen Volumes
Vcw=VA+Vc;
Vab=Vcw-Vrc-RV; 

%Alveolar Compliance
CA = VC*k*exp(-Pel*k)*((gamma + exp(-plopn/plran)*(gamma - beta*(exp(plopn/plran) + 1)))/(exp(-(Pel - plopn)/plran) + 1) -...
    exp(-plopn/plran)*(gamma - beta*(exp(plopn/plran) + 1))) +...
    (exp(-(Pel - plopn)/plran)*(gamma + exp(-plopn/plran)*(gamma - beta*(exp(plopn/plran) + 1)))*(VC - VC*exp(-Pel*k)))/(plran*(exp(-(Pel - plopn)/plran) + 1)^2);

%Resistances
Rc =  Kc*(Vcmax/Vc)^2; %Collapsible Airway Resistance
Ru = Rum+ Ku*abs(Vdot); %Upper Airway Resistance
Rs =  Rsd*exp(Ks*(VA-RV)/(TLC-RV))+Rsm; %Small Airway Resistance

prcsl=22.6974/log((exp(0.008/vrcstr)-1)/(exp(0.0001/vrcstr)-1))*prcMult; 
prctra=-prcsl*log(exp(0.008/vrcstr)-1);
pabsl=5.6/log((exp(0.0048/vabstr)-1)/(exp(0.00116/vabstr)-1)); %2.92
pabtra=-pabsl*log(exp(0.00116/vabstr)-1);

% RibCage and Abdomen Pressures
Pab=pabtra+pabsl*log(exp((Vab)/vabstr)-1);
Prc=prctra+prcsl*log(exp((Vrc)/vrcstr)-1);

Vrcdot=((Pdi_cr+Vdot*Rab+Pab-Pic-Prc)/Rrc)/(1+Rab/Rrc);

Ppl=Vrcdot*Rrc+Prc+Pdi_co+Pic;
PA=Pldyn+Ppl;
Pc=Ptm+Ppl;
Pu=Vdot*Rc+Pc;

%Airway Flows
VAdot=(Pc-PA)/Rs;
Vcdot=Vdot-VAdot;

dpdt = [(Pao-Pu-Ru*Vdot)/I; %Vdot
         VAdot/CA; %Pel
         Vcdot; %Vc
        (VAdot-Pve/Rve)/Cve; %Pve
         Vrcdot]; %Vrc
end
