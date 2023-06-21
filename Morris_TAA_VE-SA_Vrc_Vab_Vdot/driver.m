% Sensitivity analysis of TAA model, built off of Reworked Sensitivity
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set
%pars = log(pars) in DriverBasic_sens.
clear all
close all

%Set global variables for simulations
global ALLPARS INDMAP sim stateCw
Names={'TLC','RV','p_{l,opn}','p_{l,ran}','\beta',...
    '\gamma','k','v_{rc,str}','v_{ab,str}','p_{rc,sl}',...
    'RR','f','T','p_{frac}','A_{pia}',...
    'A_{mus}','R_{sm}','R_{sd}','K_s','R_{um}',...
    'K_u','I_u','p_{c,max}','p_{c,ran}','K_c',...
    'v_{c,max}','C_{ve}','R_{ve}','R_{rc}','R_{ab}'};

load Inf_Pars_TAA.mat
load Inf_Pars_VE.mat
INDMAP=union(idxparVE,idxparTAA);

%Number of breathing periods
NP=5;
%Setting which clinical outcome we are simulating
sim=9;

%Set Chest Wall compliance (low vs. high)
stateCw='low';

if sim==1 || sim==7 || sim==9
    stateCw='high';
end

%Call parameters (x0 is all parameters, allpar_names is array)
[x0,allpar_names,Init]=load_pars;


%Set identifiable parameter set
ALLPARS = x0;
pars   = x0(INDMAP); %Parameters which will undergo local sensitivity analysis
par_names=allpar_names(INDMAP,:);


[sens, yout, sols]=senseq(pars,Init,NP);  % sens come from scaled residuals
%Finds where the residuals (the data that are being used for sensitivity)
%are 0 and redefines them to be 2^-52
idx=find(yout==0);
yout(idx)=eps;
%Does the same for the parameter vector
idp=find(pars==0);
pars(idp)=eps;

% Make array for relative sensitivities
srels   = zeros(length(yout), length(pars));

%this loop weights the sensitivities
for i = 1:length(pars)
    %srels are sensitivities multiplied by parameters
    srels(:,i) = sens(:,i).*pars(i);
end

% ranked normalized relative sensitivities
[M,N] = size(sens);
sens_norm=zeros(1,N);
%Finds the l2-norm of the sensitivities for each parameter
for i = 1:N
    sens_norm(i)=norm(srels(:,i),2);
end
sens_norm_col=sens_norm';
[Rsens,Isens] = sort(sens_norm,'descend');
sens_index=([Rsens' Isens']);
sens_par_names=par_names(Isens,:);

f=figure;
plot(Rsens,'*')
set(gca,'XTick',1:numel(Rsens),'XTickLabel',sens_par_names)
set(gca,'FontSize',12);
ylabel('Relative Sensitivities')
xtickangle(60)
simstr=num2str(sim);
title(sprintf('Rsens%d',sim))
f.Position=[1 1 500 500];
saveas(gcf,sprintf('Rsens%d.fig',sim))

%return

%% Subset selection
singvals=svd(srels,0);
numopt = find(singvals./singvals(1) >= 1e-3, 1, 'last');
[cols] = gu_srrqr(srels, numopt);
    Pulm_goodpars = INDMAP(cols(1:numopt));
     Pulm_goodpars = sort([Pulm_goodpars])';
     good_par_names=allpar_names(Pulm_goodpars,:)

%savefilename='test.mat';
savefilename= ['Sens_' simstr '.mat'];
save(savefilename,'sens','par_names','Rsens','Isens','sens_par_names','sens_norm_col','singvals','numopt','Pulm_goodpars','good_par_names')