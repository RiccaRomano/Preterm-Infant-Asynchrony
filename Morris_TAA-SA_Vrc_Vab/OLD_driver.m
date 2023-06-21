%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set
%pars = log(pars) in DriverBasic_sens.
clear all
close all

global ALLPARS INDMAP
global dwMult stateCw

NP=6;

% stateCw='low';
stateCw='high';

if strcmp(stateCw,'low')==1
    dwMult=1;
else
    dwMult=0.2;
end

[x0,Init]=load_pars;
fixedpars=[8, 9, 10, 16, 17, 18]%, 28];
% these are Vstar,RV,TLC,RR,f,T,Pao
% fixedpars=[];
INDMAP  = sort(setdiff(1:length(x0), fixedpars));
ALLPARS = x0;
pars   = x0(INDMAP); %These are the ones for which we find sensitivities

[sens, yout, sols]=senseq(pars,Init,NP);  % sens come from scaled residuals
idx=find(yout==0);
yout(idx)=eps;
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
for i = 1:N
    sens_norm(i)=norm(srels(:,i),2);
end
sens_norm_col=sens_norm'
[Rsens,Isens] = sort(sens_norm,'descend');
sens_index=([Rsens' Isens']);

% savefilename= ['Sens_' stateCw 'Cw_' stateRu 'Ru_' statePao 'Pao.mat'];
savefilename=('test.mat');
save(savefilename,'sens','Rsens','Isens')
return

%% Subset selection
singvals=svd(srels,0);
numopt = find(singvals./singvals(1) >= 1e-3, 1, 'last');
[cols] = gu_srrqr(srels, numopt);
    Pulm_goodpars = INDMAP(cols(1:numopt));
     Pulm_goodpars = sort([Pulm_goodpars])';

% savefilename= [stateCw 'Cw_' stateRu 'Ru_' statePao 'Pao.mat'];
% savefilename=['test.mat']'
save(savefilename,'sens_norm_col','singvals','numopt','Pulm_goodpars')

% for i=1:size(srels,2)
% figure; plot(srels(:,i))
% end


