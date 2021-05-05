% Function to compute the POLICY FUNCTIONS and THEIR DERIVATIVES
%
%
% Inputs: - par: parameters of the policy functions
%         - fspace: structure indicating the approximation tupe, the order of approximation, etc.
%         - Grid: points on the grid of state variables where to evaluate the policy functions
%
% Outputs: All the variables of the model


function [FUNCS, States] = ComputePolFunc(par,fspace,Grid)
global  nStates nendog 
par=reshape(par,length(par)/nendog,nendog);

States     = Grid(:,1:nStates);
FUNCS      = funeval(par,fspace,Grid);

