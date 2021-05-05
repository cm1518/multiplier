%%   FUNCTION ON WHICH WE DO COLLOCATION (EULER residuals)
function res= model_colloc(par,Grid,fspace,ExoBigGrid,ExoProbMat)
global sig prob_eu
global nshocks nQuadr

nGrid = size(Grid,1);

% variables in period t
[Y,pi,Unemp,Ce,x,resid_f] = GetVar(par,fspace,Grid);

% variables in period t+1
UnempGrid = BigGrid(nQuadr^nshocks,Unemp); %expand vector, so that it can be combined with all possible realizations of future shocks
GridNext = [UnempGrid ExoBigGrid];
[Ynext,pinext,Unempnext,Cenext,xnext,resid_fnext,Gnext,anext,znext,Cunext] = GetVar(par,fspace,GridNext);

% calculate residual
RHS = (exp(-sig*anext).*((1-prob_eu).*Cenext.^(-sig) + prob_eu*Cunext.^(-sig)))./exp(pinext);
ExpRHS  = ExoProbMat * reshape(RHS,nQuadr^nshocks,nGrid); % Weight each observation by its probability

res_x    =  - x + log(ExpRHS(:));
res     = [res_x resid_f];
res     = res(:);
