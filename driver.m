%% Thoracoabdominal Asynchrony in Preterm Infant
% Author: Richard Foster and Laura Ellwein Fix
% Last edited: 3/20/2023
clear all
close all
global sim

%% Setting which clinical outcome we are simulating
% Baseline is sim=0
% High Chest Wall Compliance is sim=1
% High Upper Airway Resistance is sim=2
% High Small Airway Resistance is sim=3
% Inactive Intercostals is sim=4
% Weak Costal Diaphragm is sim=5
% Breathing Frequency is sim=6
% Combined Preterm Infant is sim=7
% Injured Lung with Low Compliance is sim=8
% Injured Lung with High Compliance is sim=9
% Obstructed Airway is sim=10
% No CPAP is sim=11
% ETT CPAP is sim=12
% Nasal CPAP is sim=13
% Nasal CPAP and Stiff Chest Wall is sim=14
% Nasal CPAP, Stiff Chest Wall, and active IA is sim=15
% Nasal CPAP, Stiff Chest Wall, active IA, and Costal Diaphragm is sim=16
sim=16;

%Number of breathing periods
NP=5;
if sim==10
    NP=10;
end

%Call parameters and initial conditions
[pars,Init]=load_pars;

Model_solver(pars,Init,NP);