function Model_solver(pars,Init,NP)
global ODE_TOLERANCE sim
load PVcurves.mat

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

VC=TLC-RV;
alpha=((1+exp(plopn/plran))*beta-gamma)/exp(plopn/plran);

tstep=0.001;
tprev=0;
times=0;
Pmussave=0;
Pdi_co_save=0;
Pdi_cr_save=0;
Picsave=0;
sols=Init;
Cdynsave=[];
Cwdynsave=[];
VTsave=zeros(NP,1);
TAAsave=zeros(NP,1);
VEsave=zeros(NP,1);
options = odeset('RelTol',ODE_TOLERANCE,'AbsTol',ODE_TOLERANCE);

for i=1:NP
    try
        if sim==10
        [Rum,Apic,f,T]=override_pars(sim,i);
        pars(12)=f;
        pars(13)=T;
        pars(15)=Apic;
        pars(20)=Rum;
        end
        tnext=tprev+T;
        tspan=tprev:tstep:tnext;
        [time,sol]=ode15s(@Model,tspan,Init,options,T,tprev,pars); %,options
        times=[times;time(2:end)];
        sols=[sols; sol(2:end,:)];
        Init=[sol(end,1:5)];
        Peltemp=sol(2:end,2);
        maxPel=max(Peltemp);
        minPel=min(Peltemp);
        endPel=Peltemp(end);
        Frectemp=alpha+(gamma-alpha)./(1+exp(-(Peltemp-plopn)/plran));
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
        Vrctemp=sol(2:end,5);
        Vabtemp=Vcwtemp-Vrctemp-RV;
        imaxPrc=find(abs(Vcw_range-Vcwmax)<0.0001, 1, 'last' );
        iminPrc=find(abs(Vcw_range-Vcwmin)<0.0001, 1, 'last' );
        maxPrc=Pel_range(imaxPrc);
        minPrc=Pel_range(iminPrc);
        Cwdyn=(Vcwmax-Vcwmin)/(maxPrc-minPrc);
        Cwdynsave=[Cwdynsave;Cwdyn];
        VT=VAmax-VAmin;
        VTsave(i)=VT*1000;
        VE=VT*f*60;
        VEsave(i)=VTsave(i)*f*60;
        TAAsave(i)=pk(Vrctemp,Vabtemp,T,tstep);
        Pmustemp=(Amus*cos(2*pi*f*(tspan-tprev)) + -Amus)';
        Pmussave=[Pmussave; Pmustemp(2:end)];
        Pdi_co_temp=Pfrac*Pmustemp;
        Pdi_co_save=[Pdi_co_save; Pdi_co_temp(2:end)];
        Pdi_cr_temp=Pmustemp-Pdi_co_temp;
        Pdi_cr_save=[Pdi_cr_save; Pdi_cr_temp(2:end)];
        Pictemp=(Apic*cos(2*pi*f*(tspan-tprev)) + -Apic)';
        Picsave=[Picsave; Pictemp(2:end)];
        tprev = tnext;
        disp([i,VAmax*1000,VAmin*1000,VT*1000,f*60,VE*1000,Cdyn*1000,Cwdyn*1000,fracOpen,minPel])
    catch
        disp([i,VAmax*1000,VAmin*1000,VT*1000,f*60,VE*1000,Cdyn*1000,Cwdyn*1000,fracOpen,minPel]) %Cwdyn*1000,
        time=[]; sol=[];
        break
    end
end

Vdot=sols(:,1);  Pel=sols(:,2);  Vc=sols(:,3); Pve=sols(:,4); Vrc=sols(:,5);

Frec=alpha+(gamma-alpha)./(1+exp(-(Pel-plopn)./plran));
Vel=VC*(1-exp(-k*Pel));
VA=Frec.*Vel+RV;
Vcw=VA+Vc;
Vab=Vcw-Vrc-RV;
Ru = transpose(Rum+ Ku*abs(Vdot));
Rc =  Kc*(Vcmax./Vc).^2;
Rs = Rsd*exp(Ks*(VA-RV)./(TLC-RV))+Rsm;
Rtot=Rc+Ru+Rs;
prcsl=22.6974/log((exp(0.008/vrcstr)-1)/(exp(0.0001/vrcstr)-1))*prcMult; 
prctra=-prcsl*log(exp(0.008/vrcstr)-1);
pabsl=5.6/log((exp(0.0048/vabstr)-1)/(exp(0.00116/vabstr)-1))*1; %2.92
pabtra=-pabsl*log(exp(0.00116/vabstr)-1);

Pab=pabtra+pabsl*log(exp(Vab/vabstr)-1);
Prc=prctra+prcsl*log(exp(Vrc/vrcstr)-1);
Pldyn=Pel+Pve;
Ptm = pcmax-pcran*log(Vcmax./(Vc-0)-1);

Vrcdot=((Pdi_cr_save+Vdot*Rab+Pab-Picsave-Prc)/Rrc)/(1+Rab/Rrc);

Ppl=Vrcdot*Rrc+Prc+Pdi_co_save+Picsave;
Pc=Ptm+Ppl;
PA=Pldyn+Ppl;
Pu=Vdot.*Rc+Pc;

t=times;

simstr=num2str(sim);
savefilename= ['Sim' simstr '.mat'];
save(savefilename);

disp('Phase Angle is: ');
disp(TAAsave(end));
disp('Tidal Volume is: ');
disp(VTsave(end));
 
%% Figures
% start=101;
% figure(1)
% hold on
% idx=find(VA_range==min(VA_range));
% plot(Pel_range,Vcw_range*1000,'-','Color',[0 0 0],'LineWidth',2); %Chest Wall Compliance
% plot(Pel_range(idx:end),VA_range(idx:end)*1000,'-','Color',[0.75 0.75 0.75],'LineWidth',2); %Lung Compliance
% plot(Ptot,VolN*1000,'-.','Color',[0.5 0.5 0.5],'LineWidth',2); %Total Repiratory Compliance
% plot(Pldyn(start:end),VA(start:end)*1000,'-','Color',[0 0 0],'LineWidth',2) %Tidal Loop
% plot(P_FRC,FRC*1000,'o','Color',[0 0 0],'LineWidth',2);
% plot(-P_FRC,FRC*1000,'o','Color',[0 0 0],'LineWidth',2); %Points of Intersection with FRC
% line([-5 35],[FRC*1000 FRC*1000],'Color',[0.5 0.5 0.5],'Linestyle',':','LineWidth',2);
% line([0 0],[0 60],'Color',[.7 .7 .7]);
% set(gca,'FontSize',12);
% title('Simulation 0','FontSize',18);
% legend('Chest Wall','Lung','Respiratory','Tidal loop','','','FRC','','Location','Southeast');
% axis([-5 35 RV*1000 60]);
% 
% figure(2)
% plot(t,Ru,t,Rs,t,Rc);
% title('Resistance');
% legend('Upper Airways','Small Airways','Collapsible Airways');
% 
% figure(5)
% plot(t,Vrc*1000,'-b',t,Vab*1000,'-r');
% legend('Rib Cage','Abdomen');
% 
% figure(6)
% plot(t,Vrcdot,t,Vdot-Vrcdot)
% legend('Vrcdot','Vabdot')
% 
% figure(7)
% plot(t,Prc,t,Pab,t,Ppl,t,Pmussave);
% legend('Prc','Pab','Ppl','Pmus');