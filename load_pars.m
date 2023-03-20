function [pars,Init]=load_pars
global ODE_TOLERANCE DIFF_INC sim

ODE_TOLERANCE=1e-8;
DIFF_INC = sqrt(ODE_TOLERANCE);

TLC   = .001*63; %Total Lung Capacity
RV    = .001*23; %Residual Volume
plopn =0.1; %p_l,opn, Mean Opening Lung Pressure
plran = 0.5; %p_l,ran, Lung pressure range coefficient
beta=0.01; %lung Recruitment baseline fraction
gamma = 1; %Lung Recruitment Maximum Fraction
k=0.05; %Pressure-Dependence Lung Coefficient
vrcstr=0.00425; %RC PV-Curve Stretch Coefficient
vabstr=0.004; %AB PV-Curve Streth Coefficient
prcMult=0.4; %RC PV-Curve Slope Multiplier
RR      = 60; %Respiratory Rate
f       = RR/60;  %Frequency
T       = 1/f; %Period
Pfrac =.52; %portion of Pmus going to Costal Diaphragm
Apic=1; %Intercostal/Accessory Amplitude
Amus=4.4; %Diaphragm Muscle Amplitude
Rsm    = 12; %Small Airway Resistance Minimum
Rsd    = 20; %Small Airway Resistance Difference
Ks    = -15; %Small Airway Volume-Dependence Coefficient
Rum   = 20; %Upper Airway Resistance Minimum Value
Ku    = 60; %Upper Airway Resistance Flow-Dependence Coefficient
I     = 0.33; %Upper Airway Inertance
pcmax    = 4.4; %Collapsible Airway PV-Curve Pressure at Max Compliance
pcran   = 4.4; %Collapsible Airway PV-Curve Pressure Range Coefficient
Kc    = .1; %Collapsible Airway Volume-Dependence Coefficient
Vcmax = 0.0025; %Collapsible Airway Max Volume
Cve   = 0.005; %Viscoelastic Lung Compliance
Rve   = 20; %Viscoelastic Lung Resistance
Rrc   = 1; %Rib Cage Resistance
Rab   = 1; %Abdomen Resistance
Pao=0; %Airway Opening Pressure


switch sim
    case 0
        disp('baseline')
    case 1
        disp('high Cw')
        prcMult=0.08;
    case 2
        disp('high Rum')
        Rum=230;
    case 3
        disp('high Rsm')
        Rsm=120;
    case 4
        disp('0 Apic')
        Apic=0;
    case 5
        disp('low Pfrac')
        Pfrac=0.28;
    case 6
        disp('high RR')
        RR=100;
        f=RR/60;
        T=1/f;
    case 7
        disp('Preterm')
        prcMult=0.08;
        Rum=230;
        Rsm=120;
        Apic=0;
        Pfrac=0.28;
        RR=100;
        f=RR/60;
        T=1/f;
    case 8
        disp('Injured lung, low Cw')
        prcMult=0.4;
        plopn=15;
        plran=3;
        beta=0.1;
    case 9
        disp('Injured lung, high Cw')
        prcMult=0.08;
        plopn=15;
        plran=3;
        beta=0.1;
    case 10
        disp('Obstructed airway')
    case 11
        disp('Injured Lung with Combined Infant, No CPAP')
        prcMult=0.08;
        plopn=15;
        plran=3;
        beta=0.1;
        gamma=1;
        Rsm=120;
        Apic=0;
        Pfrac=0.28;
        RR=60;
        f=RR/60;
        T=1/f;
        Pao=0;
    case 12
        disp('ETT-CPAP');
        prcMult=0.08;
        Rum=100;
        Rsm=120;
        Apic=0;
        Pfrac=0.28;
        RR=60;
        f=RR/60;
        T=1/f;
        plopn=15;
        plran=3;
        beta=0.1;
        gamma=1;
        Pao=8;
    case 13
        disp('Nasal-CPAP')
        prcMult=0.08;
        plopn=15;
        plran=3;
        beta=0.1;
        gamma=1;
        Rsm=120;
        Apic=0;
        Pfrac=0.28;
        RR=60;
        f=RR/60;
        T=1/f;        
        Pao=8;
    case 14
        disp('Nasal-CPAP and Stiff chest wall');
        prcMult=0.4;
        plopn=15;
        plran=3;
        beta=0.1;
        gamma=1;
        Rsm=120;
        Apic=0;
        Pfrac=0.28;
        RR=60;
        f=RR/60;
        T=1/f;        
        Pao=8;
    case 15
        disp('Nasal-CPAP and Stiff chest wall, active IA');
        prcMult=0.4;
        plopn=15;
        plran=3;
        beta=0.1;
        gamma=1;
        Rsm=120;
        Apic=1;
        Pfrac=0.28;
        RR=60;
        f=RR/60;
        T=1/f;        
        Pao=8;
    case 16
        disp('Nasal-CPAP and Stiff chest wall, active IA, and strong costals');
        prcMult=0.4;
        plopn=15;
        plran=3;
        beta=0.1;
        gamma=1;
        Rsm=120;
        Apic=1;
        Pfrac=0.52;
        RR=60;
        f=RR/60;
        T=1/f;        
        Pao=8;
    otherwise
        disp('inappropriate output')
end

curvepars=[TLC RV plopn plran beta gamma k vrcstr vabstr prcMult];
P_FRC=test_curves_preterm_freeRCAB(curvepars);
Pel0=P_FRC;
Vdot0=0;
Vc0 = 0.0001;
Pve0 = 0;
Vrc0=0.0001;
Init=[Vdot0 Pel0 Vc0 Pve0 Vrc0]; 

pars=[TLC RV plopn plran beta gamma k vrcstr vabstr prcMult RR f T Pfrac Apic Amus Rsm Rsd Ks Rum Ku I pcmax pcran Kc Vcmax Cve Rve Rrc Rab Pao]';
