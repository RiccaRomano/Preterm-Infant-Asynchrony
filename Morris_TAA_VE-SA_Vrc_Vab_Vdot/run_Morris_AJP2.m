%% Morris screening
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES

clc; clear; 
global sim

sim = 0

Names={'TLC','RV','p_{l,max}','p_{l,ran}','\beta',...
    '\gamma','k','v_{rc,str}','v_{ab,str}','p_{rc,str}',...
    'RR','f','T','p_{frac}','A_{pia}',...
    'A_{mus}','R_{sm}','R_{sd}','K_s','R_{um}',...
    'K_u','I_u','p_{c,max}','p_{c,ran}','K_c',...
    'v_{c,max}','C_{ve}','R_{ve}','R_{rc}','R_{ab}'};


%% Initialize things for the model
NC     = 60;                        % Number of Morris samples
[pars,~,Init]=load_pars;

dt     = .002;
T = 1;
tspace = 0:dt:T;

fixedpars=[1,2,11,12,13];
par_ids  = sort(setdiff(1:length(pars), fixedpars));
par_char=Names(par_ids);
par_nom=pars;

% Define a function call that gives you the outputs (pressure,flow, area)
f = @(q) Model_solver_Morris(q,Init,tspace,par_ids,par_nom);
%% Morris screening parameters
q_nom = pars(par_ids); % Nominal parameter to perturb around
r = 120;               % Max samples to try (e.g., r-smallR is how many failures we can account for)
smallR = 100;          % Number of samples we want to actually use (smaller than r)
p = length(par_ids);   % Number of parameters
l = 60;                % Number of levels
delta = l./(2*(l-1));  % Step Size
% MJC 9/13/2021 - new file to get more parameter specific bounds:
[upper,lower] = get_bounds_AJP(par_nom,par_ids);

d = cell(r,p);
%% Try to use the randomization algorithm
% Note that all parameters are scaled to be in the range 0,1 and then
% rescaled in the model evaluation.
A = zeros(p+1,p);
for i=1:p
    A(i+1:p+1,i) = ones((p-(i-1)),1);
end
X = zeros((p+1).*r,p.*r);
%% Begin randomization algorithm, see Smith UQ book for details
%Trajectory Design or Winding Stairs method for sampling points
F_storage = cell(p+1,smallR);
qstar = unifrnd(0,1-delta,r,p);
Jp = ones(p+1,p);
J1 = ones(p+1,1);
P = eye(p,p);
UL_MAT = eye(p).*(upper-lower);
i=1;
func_evals = 1;
while i<=smallR
    i
    exitflag = 0;
    if func_evals == r
        error('Parameter Space Infeasible: Exiting.\n');
    end
    qcurr = qstar(i,:);
    pm1 =rand(p,1);
    Dstar = eye(p).*(pm1 > 0.5) - eye(p).*(pm1 <= 0.5);
    [who,where] = sort(rand(p,1));
    Pstar = P(where,:);
    Astar = (J1*qcurr + (delta./2).*((2.*A - Jp)*Dstar + Jp))*Pstar;
    C = Astar*UL_MAT+J1*lower;
    fpast = f(C(1,:));
    F_storage{1,i} = fpast;
    for j=1:p
        disp([j+1 Names(par_ids(where(j)))])
        fpresent = f(C(j+1,:));
        if isempty(fpresent) || isempty(fpast)
           for s=j:-1:1
              d{i,where(s)} = {}; %Clear all previous entries 
           end
           exitflag = 1;
           break;
        end
        d{i,where(j)} = (fpresent - fpast)./delta;

        fpast = fpresent;
        F_storage{where(j)+1,i} = fpast;
    end
    if exitflag == 1
        exitflag = 0;
    else
        i = i+1;
    end
end
save('Morris_results_TAA_VE');