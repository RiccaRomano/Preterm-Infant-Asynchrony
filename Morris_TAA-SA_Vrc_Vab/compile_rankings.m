% I will load in each file Sens_#.mat which contains the variables sens,
% Rsens, Isens, and sens_par_names.
% sens is the sensitivity matrix with states scaled by the max value of the
% state (not multiplied by parameter value). Rsens is ranked normed sensitivity scalars, and Isens is the
% associated index of the parameter vector. 

% Just realized this whole time rout was only 1 of the two states.

% I want to create an array of unranked sensivities. Each row is a
% parameter progressing as in the order of the parameter vector. Each
% column is a simulation. Then I can create some sort of average ranking
% and average sensitivity, but display as row numbers on a chart.
% Whichever parameters are always insensitive can be fixed in subset
% selection and covariance analyses, and potentially be used to reduce the model.

% Then do Morris elementary effects with all parameters on TAA and tidal
% volume. 