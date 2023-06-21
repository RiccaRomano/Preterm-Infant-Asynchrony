
%All volume inputs are in absolute units, L
global Rum Ku Kc Vcmax Rsd Ks Rsm Vstar RV TLC 
global Cve Rve  I cc dc VC  %volume0Pres
global RR f T Phalf Ptau k  
global beta gamma alpha 
global dwMult Rrc Rab 
global brc crc drc bab cab dab

Phalf =.1;%8;%
Ptau =.5;%2;%
beta=0.01;%.1;%
gamma = 1;
alpha=((1+exp(Phalf/Ptau))*beta-gamma)/exp(Phalf/Ptau);
k=0.05;%0.05;

% RV     =0.023;
% TLC   = 0.063;
VC  = TLC-RV;
% alpha_cw=0.25;%
% volume0Pres = alpha_cw*VC+RV;
RR      = 60;
f       = RR/60;  
T       = 1/f;

cc    = 4.4;
dc    = 4.4;
Rum   = 20;
Ku    = 60;

I     = 0.33;
Kc    = .1;
Vcmax = 0.0025;
Rsm    = 12;
Rsd    = 20;
Ks    = -15;
Vstar = TLC;
Cve   = 0.005*1;
Rve   = 20;
Rrc   = 1;
Rab   = 1;

%% Free Rib Cage and Abdominal Constants
brc=0.00425 %normal bw=0.00425
drc=22.6974/log((exp(0.008/brc)-1)/(exp(0.0001/brc)-1))*dwMult 
crc=-drc*log(exp(0.008/brc)-1)

bab=0.004
dab=5.6/log((exp(0.0048/bab)-1)/(exp(0.00116/bab)-1))*1 %2.92
cab=-dab*log(exp(0.00116/bab)-1)


