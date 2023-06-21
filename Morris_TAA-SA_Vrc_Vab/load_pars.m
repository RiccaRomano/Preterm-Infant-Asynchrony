function [pars,par_names,Init]=load_pars
global ODE_TOLERANCE DIFF_INC
global stateCw sim

ODE_TOLERANCE=1e-8;
DIFF_INC = sqrt(ODE_TOLERANCE);

TLC   = .001*63;
RV    = .001*23;
Phalf =.1; %cF Y
Ptau = 0.5; %dF Y
beta=0.01; %Y
gamma = 1; %Y
k=0.05; %Y
brc=0.00425;  %Y
bab=0.004; %Y
if strcmp(stateCw,'low')==1
    drcMult=0.4;
else
    drcMult=0.4/5;
end

curvepars=[TLC RV Phalf Ptau beta gamma k brc bab drcMult];
P_FRC=test_curves_preterm_freeRCAB(curvepars);
Pel0=P_FRC;
Vdot0=0;
Vc0 = 0.0001;
Pve0 = 0;
Vrc0=0.0001;
Init=[Vdot0 Pel0 Vc0 Pve0 Vrc0]; 

RR      = 60; %Y
f       = RR/60;  %Y 
T       = 1/f; %Y

Pfrac =.52; %fraction of Pmus going to costal (main) diaphragm %Y
Apic=1; %Y
Amus=4.4;% %Y

Rsm    = 12; %Y
Rsd    = 20; %Y
Ks    = -15; %Y
Rum   = 20; %Y
Ku    = 60; %Y
I     = 0.33; %Y
cc    = 4.4; %Y
dc    = 4.4; %Y
Kc    = .1; %Y
Vcmax = 0.0025; %Y
Cve   = 0.005; %Y
Rve   = 20; %Y
Rrc   = 1; %Y
Rab   = 1; %Y

switch sim
    case 0
        disp('baseline')
    case 1
        disp('high Cw')
    case 2
        disp('high Rum')
        Rum=230;
    case 3
        disp('high Rsm')
        Rsm=120;
    case 4
        disp('0 Apic')
        Apic=0.01;
    case 5
        disp('low Pfrac')
        Pfrac=0.28;
    case 6
        disp('high RR')
        RR=100;
    case 7
        disp('Preterm')
        Rum=230;
        Rsm=120;
        Apic=0.01;
        Pfrac=0.28;
        RR=100;
    case 8
        disp('Unhealthy lung, low Cw')
        Phalf=8;
        Ptau=3;
        beta=0.1;
    case 9
        disp('Unhealthy lung, high Cw')
        Phalf=8;
        Ptau=3;
        beta=0.1;
    case 10
        disp('Obstructed airway')
        Rum=2000;
    otherwise
        disp('inappropriate output')
end

pars=[TLC RV Phalf Ptau beta gamma k brc bab drcMult RR f T Pfrac Apic Amus Rsm Rsd Ks Rum Ku I cc dc Kc Vcmax Cve Rve Rrc Rab]';
par_names={'TLC','RV','p_{l,opn}','p_{l,ran}','\beta',...
    '\gamma','k','v_{rc,str}','v_{ab,str}','p_{rc,sl}',...
    'RR','f','T','p_{frac}','A_{pia}',...
    'A_{mus}','R_{sm}','R_{sd}','K_s','R_{um}',...
    'K_u','I_u','p_{c,max}','p_{c,ran}','K_c',...
    'v_{c,max}','C_{ve}','R_{ve}','R_{rc}','R_{ab}'}';

%(1,2,11,12,13) are TLC,RV,RR,f,T
%l=length(pars);
