function [upper,lower] = get_bounds_AJP(pars,ids)
% pars to screen: 3:10, 14:30
npar = length(pars);
lower = ones(1,npar);
upper = ones(1,npar);

lmult=.75;
umult=1.25;

lower=pars(ids)'*lmult;
upper=pars(ids)'*umult;

%% MY BOUNDS
% lower(3)=0.1;   upper(3)=0.5; %Phalf
% lower(4)=0.5;     upper(4)=5; %Ptau
% lower(5)=.1;   upper(5)=1; %beta
% lower(6)=.6;    upper(6)=1; %gamma
% lower(7)=.04;   upper(7)=0.1; %k
% lower(8)=0.001; upper(8)=0.01; %brc
% lower(9)=0.001; upper(9)=0.01; %bab
% lower(10)=.04;  upper(10)=.8; %drcMult
% lower(14)=.2;   upper(14)=.8; %Pfrac
% lower(15)=0.001;    upper(15)=1; %Apic
% lower(16)=1;    upper(16)=5; %Amus
% lower(17)=10;  upper(17)=20; %Rsm
% lower(18)=10;  upper(18)=30; %Rsd
% lower(19)=-20;  upper(19)=-5; %Ks
% lower(20)=10;   upper(20)=30; %Rum
% lower(21)=40; upper(21)=70; %Ku
% lower(22)=0.1;  upper(22)=0.5; %I
% lower(23)=3;   upper(23)=5; %cc
% lower(24)=3;   upper(24)=5;%dc
% lower(25)=0.1;   upper(25)=0.5; %Kc
% lower(26)=.001; upper(26)=.004; %Vcmax
% lower(27)=.001; upper(27)=0.01; %Cve
% lower(28)=10;   upper(28)=40; %Rve
% lower(29)=0.1;  upper(29)=10; %Rrc
% lower(30)=0.1;  upper(30)=10; %Rab



%% DRIVE BOUNDS
lower(3)=0.1;   upper(3)=15; %Phalf GOOD
lower(4)=0.5;     upper(4)=3; %Ptau GOOD
lower(5)=0.01;   upper(5)=1; %beta GOOD
lower(6)=0.6;    upper(6)=1; %gamma GOOD
lower(7)=0.04;   upper(7)=0.07; %k GOOD
lower(8)=0.001; upper(8)=0.01; %brc GOOD
lower(9)=0.001; upper(9)=0.01; %bab GOOD
lower(10)=0.04;  upper(10)=0.8; %drcMult GOOD
lower(14)=0.2;   upper(14)=0.8; %Pfrac GOOD
lower(15)=0.001;    upper(15)=1; %Apic GOOD
lower(16)=1;    upper(16)=5; %Amus GOOD
lower(17)=10;  upper(17)=120; %Rsm 
lower(18)=10;  upper(18)=30; %Rsd
lower(19)=-20;  upper(19)=-10; %Ks GOOD
lower(20)=10;   upper(20)=300; %Rum GOOD
lower(21)=40; upper(21)=70; %Ku
lower(22)=0.1;  upper(22)=0.5; %I
lower(23)=3;   upper(23)=5; %cc GOOD
lower(24)=3;   upper(24)=5;%dc GOOD
lower(25)=0.1;   upper(25)=10; %Kc
lower(26)=.002; upper(26)=.003; %Vcmax
lower(27)=.001; upper(27)=0.01; %Cve
lower(28)=10;   upper(28)=40; %Rve
lower(29)=0.1;  upper(29)=10; %Rrc GOOD
lower(30)=0.1;  upper(30)=10; %Rab GOOD

upper = upper(ids);%.*pars(ids);
lower = lower(ids);%.*pars(ids);