
% increase the size of the variables to account for all the possibile
% realizations of the shock.


function varargout = BigGrid(S,varargin)
ninputs = nargin-1; 
for i=1: ninputs
    varargout{i} = kron(varargin{i},ones(S,1));
end