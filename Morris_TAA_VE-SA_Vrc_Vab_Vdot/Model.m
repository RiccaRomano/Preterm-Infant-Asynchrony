function [dpdt] =Model(t,p,T,tprev,pars)

TLC=pars(1);
RV=pars(2);
Phalf=pars(3);
Ptau=pars(4);
beta=pars(5);
gamma=pars(6);
k=pars(7);
brc=pars(8);
bab=pars(9);
drcMult=pars(10);
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
cc=pars(23);
dc=pars(24);
Kc=pars(25);
Vcmax=pars(26);
Cve=pars(27);
Rve=pars(28);
Rrc=pars(29);
Rab=pars(30);

Vdot = p(1);
Pel = p(2);
Vc = p(3);
Pve = p(4);
Vrc = p(5);

Pldyn=Pel+Pve;
Ptm = cc-dc*log(Vcmax./(Vc-0)-1);

% Driving Pressures
%%% Toggle for muscle pressure drive function
Pmus = (Amus*cos(2*pi*f*(t-tprev)) + -Amus); %
Pdi_co=Pfrac*Pmus; %Pfrac is defined in the parameters script
Pdi_cr=Pmus-Pdi_co;
Pic = (Apic*cos(2*pi*f*(t-tprev)) + -Apic);

VC  = TLC-RV;
alpha=((1+exp(Phalf/Ptau))*beta-gamma)/exp(Phalf/Ptau);
Frec=alpha+(gamma-alpha)./(1+exp(-(Pel-Phalf)/Ptau));
Vel=VC*(1-exp(-k*Pel));
VA=Frec*Vel+RV;

Vcw=VA+Vc;
Vab=Vcw-Vrc-RV; 

CA = VC*k*exp(-Pel*k)*((gamma + exp(-Phalf/Ptau)*(gamma - beta*(exp(Phalf/Ptau) + 1)))/(exp(-(Pel - Phalf)/Ptau) + 1) -...
    exp(-Phalf/Ptau)*(gamma - beta*(exp(Phalf/Ptau) + 1))) +...
    (exp(-(Pel - Phalf)/Ptau)*(gamma + exp(-Phalf/Ptau)*(gamma - beta*(exp(Phalf/Ptau) + 1)))*(VC - VC*exp(-Pel*k)))/(Ptau*(exp(-(Pel - Phalf)/Ptau) + 1)^2);

% Resistances
Rc =  Kc*(Vcmax/Vc)^2;
Ru = Rum+ Ku*abs(Vdot);
Rs =  Rsd*exp(Ks*(VA-RV)/(TLC-RV))+Rsm;

drc=22.6974/log((exp(0.008/brc)-1)/(exp(0.0001/brc)-1))*drcMult; 
crc=-drc*log(exp(0.008/brc)-1);
dab=5.6/log((exp(0.0048/bab)-1)/(exp(0.00116/bab)-1))*1; %2.92
cab=-dab*log(exp(0.00116/bab)-1);

% RibCage and Abdomen Pressures
%Pab and Prc are now functions of %VC, we need to make them in terms of
%absolute volume again
Pab=cab+dab*log(exp((Vab)/bab)-1);
Prc=crc+drc*log(exp((Vrc)/brc)-1);

Vrcdot=((Pdi_cr+Vdot*Rab+Pab-Pic-Prc)/Rrc)/(1+Rab/Rrc);

Ppl=Vrcdot*Rrc+Prc+Pdi_co+Pic;
PA=Pldyn+Ppl;
Pc=Ptm+Ppl;
Pu=Vdot*Rc+Pc;

VAdot=(Pc-PA)/Rs;
Vcdot=Vdot-VAdot;

dpdt =   [(0-Pu-Ru*Vdot)/I;              %Vdot
    (VAdot)/CA; ...                %Pel
    Vcdot;                         %Vc
    (VAdot-Pve/Rve)/Cve;          %Pve
        Vrcdot];  %Vrc

end

%%%%

% if tprev>tcpap
%    Pao=5;%5;
% else
%Pao=0;
% end
