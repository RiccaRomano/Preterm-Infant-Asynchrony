function [sols,rout]=Model_solver(pars,Init,NP)
% close all
 
global ODE_TOLERANCE
load PVcurves.mat

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

VC=TLC-RV;
alpha=((1+exp(Phalf/Ptau))*beta-gamma)/exp(Phalf/Ptau);

tstep=0.001;
tprev=0;
tnext=tprev+T;
periods=[tprev tnext];
tspan = [tprev:tstep:tnext]; %Creates timespan over which to solve DEs
ppp=length(tspan);


ts=0;
times=0;
sols=Init;

Pmussave=0;
Pdi_co_save=0;
Pdi_cr_save=0;
Picsave=0;
VTsave=0.005;
VAminsave=[];
Cdynsave=[];
Cwdynsave=[];

fracAtel=0;

% options = [];%odeset('OutputFcn',@outputFcn);
% Solve ODEs from tstart to TFinal
options = odeset('RelTol',ODE_TOLERANCE,'AbsTol',ODE_TOLERANCE);

for i=1:NP
    try
        [time,sol]=ode15s(@Model,tspan,Init,options,T,tprev,pars); %,options
        times=[times;time(2:end)];
        sols=[sols; sol(2:end,:)];
        Init=[sol(end,1:5)];

        Peltemp=sol(2:end,2);
        maxPel=max(Peltemp);
        minPel=min(Peltemp);
        endPel=Peltemp(end);

        Frectemp=alpha+(gamma-alpha)./(1+exp(-(Peltemp-Phalf)/Ptau));
        fracOpen=max(Frectemp);
        minFrec=min(Frectemp);
        endFrec=Frectemp(end);
        Veltemp=VC*(1-exp(-k*Peltemp));
        VAtemp=Frectemp.*Veltemp+RV;
        VAmax=max(VAtemp);
        VAmin=min(VAtemp);

        Cdyn=(VAmax-VAmin)/(maxPel-minPel);
        Cdynsave=[Cdynsave;Cdyn];
        Vcwtemp=VAtemp+sol(2:end,3);
        Vcwmax=max(Vcwtemp);
        Vcwmin=min(Vcwtemp);
        
        imaxPrc=find(abs(Vcw_range-Vcwmax)<0.0001, 1, 'last' );
        iminPrc=find(abs(Vcw_range-Vcwmin)<0.0001, 1, 'last' );
        maxPrc=Pel_range(imaxPrc);
        minPrc=Pel_range(iminPrc);
        Cwdyn=(Vcwmax-Vcwmin)/(maxPrc-minPrc);
        Cwdynsave=[Cwdynsave;Cwdyn];

        VT=VAmax-VAmin;
        VTsave=[VTsave;VT];
        VE=VT*f*60;

        Pmustemp=(Amus*cos(2*pi*f*(tspan-tprev)) + -Amus)';
        Pmussave=[Pmussave; Pmustemp(2:end)];

        tprev = tnext;
        tnext=tprev+T;
        tspan = [tprev:tstep:tnext];
        ppp=length(tspan);
        pend=periods(end)+T;
        periods=[periods pend];
    catch
        disp([i,VAmax*1000,VAmin*1000,VT*1000,f*60,VE*1000,Cdyn*1000,Cwdyn*1000,fracOpen,minPel]) %Cwdyn*1000,
        time=[]; sol=[];
        break
    end
end

t=times;
starti=round(T/tstep*1.09);
endi=round(starti+T/tstep);
Vdot=sols(:,1);  Pel=sols(:,2);  Vc=sols(:,3); Pve=sols(:,4); Vrc=sols(:,5);
 % PelR=sols(starti:endi,2);  VcR=sols(starti:endi,3); PveR=sols(starti:endi,4);

Frec=alpha+(gamma-alpha)./(1+exp(-(Pel-Phalf)./Ptau));
Vel=VC*(1-exp(-k*Pel));
VA=Frec.*Vel+RV;

Vcw=VA+Vc;
Vab=Vcw-Vrc-RV;

for k=1:length(Vdot)
%     if Vdot(k)<0
%         Ru(k) = (Rum+ Ku*abs(Vdot(k)))*10;%
%     else
        Ru(k) = Rum+ Ku*abs(Vdot(k));
%     end
end
Ru=Ru';
Rc =  Kc*(Vcmax./Vc).^2;
Rs = Rsd*exp(Ks*(VA-RV)./(TLC-RV))+Rsm;
Rtot=Rc+Ru+Rs;

drc=22.6974/log((exp(0.008/brc)-1)/(exp(0.0001/brc)-1))*drcMult; 
crc=-drc*log(exp(0.008/brc)-1);
dab=5.6/log((exp(0.0048/bab)-1)/(exp(0.00116/bab)-1))*1; %2.92
cab=-dab*log(exp(0.00116/bab)-1);

Pab=cab+dab*log(exp(Vab/bab)-1);
Prc=crc+drc*log(exp(Vrc/brc)-1);
Pldyn=Pel+Pve;
Ptm = cc-dc*log(Vcmax./(Vc-0)-1);

Vrcdot=((Pdi_cr_save+Vdot*Rab+Pab-Picsave-Prc)/Rrc)/(1+Rab/Rrc);

Ppl=Vrcdot*Rrc+Prc+Pdi_co_save+Picsave;
Pc=Ptm+Ppl;
PA=Pldyn+Ppl;
Pu=Vdot.*Rc+Pc;

% PldynR = PelR+PveR;

VdotR=sols(starti:endi,1);
% PplR=Ppl(starti:endi);
%  rout=[VdotR/max(abs(VdotR));PplR/max(abs(PplR))]; %
 VabR=Vab(starti:endi);
 VrcR=Vrc(starti:endi);
 rout=[VabR/max(abs(VabR));VrcR/max(abs(VrcR));VdotR/max(abs(VdotR))];
%rout=[VdotR/max(abs(VdotR))]; %this is an experiment
% rout=VdotR/max(abs(VdotR));
% rout=[VdotR;PplR];

