function [sols,rout] = Model_wrapper(pars,Init,NP)
% Wrapper function to account for parameters that are fixed

global INDMAP ALLPARS 

tempars = ALLPARS;
tempars(INDMAP) = pars;

[sols,rout] = Model_solver(tempars,Init,NP);


end
