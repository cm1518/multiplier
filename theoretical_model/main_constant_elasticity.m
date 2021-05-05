addpath(strcat(pwd,'/CEtools'));


%% ----
% Main code for "Underdstanding the Size of the Government Spending Multiplier: It's in the Sign)
% Barnichon, Debortoli, Matthes (ResStud, 2020)
% This version: Sept. 2020.

clear vars; clc; close all;

%set seed
rng(10);

global phipi sig gamma0 gamma1 gamma2 prob_eu Ce_ss Cu_ss G_ss U_ss
global nshocks nQuadr nStates nendog do_fullins do_constantelast

%% parameters of the model
betta   = 0.99;             % discount factor
sig     = 2;                % risk aversion
alpa    = 0;                % 1 - labor share
phiy    = 0;                % coefficient on output-gap in Taylor rule
gbar    = 0.2;              % G/Y ratio in steady state (data: 16% only G, 26% if including gov. investm)
prob_eu =  1-(1-0.03)^3;    % separation rate (based on 3% monthly separation rate)
U_ss    = 0.045;            % steady state unemployment (prob_eu/(prob_eu+prob_ue));

theta   = 0.4;              % target for net replacement rate b/(W-T) in steady state
phipi   = 1.3;              % coefficient on inflation in Taylor rule

rho_g   = 0.8;              % persistence of government expenditure shock
sigma_epsg = 0.010;         % standard deviation of government expenditure shock

rho_z   = 0.18;             % persistence of demand shock
sigma_epsz = 0.05;          % standard deviation of demand shock

%% parameters of solution method
nQuadr  = 11;               % number of quadrature points to calculate expectations 
nshocks = 2;                % number of exogenous shocks
nStates = 3;                % number of state variables
nendog  = 2;                % number of policy functions to be solved

ApproxType='cheb'; splineorder=[];  % type of polynominal approximation
Order = [2 7 7];                    % order of approximation
Umin = U_ss;                        % lower bound on grid for unemployment
Umax = U_ss + 0.05;                 % upper bound on grid for unemployment

do_fullins = 0;             % set this to 1 to solve with perfect risk-sharing
do_constantelast = 1;       % set this to 1 to solve model with constant wage elasticity



%% Calculate steady state values
G_ss    = gbar*(1-U_ss);
Ce_ss   = (1-U_ss-G_ss)/(theta*U_ss + (1-U_ss));
Cu_ss   = theta*Ce_ss;% - G_ss;

if do_constantelast
    gamma2 = 1.4419;
else
    gammas = fsolve('calib_gammas',[0.99;1;5]);
    calib_gammas(gammas);
    gamma0   = gammas(1);        % level of wage rigidity
    gamma1   = gammas(2);        % slope of wage rigidity
    gamma2   = gammas(3);        % curvature of wage rigidity
    
    unemp_plot = [U_ss:0.01:U_ss+0.03];
    plot(unemp_plot,400*log(rigidfunc(unemp_plot,U_ss)));
    elast =  gamma1^(-gamma2)*U_ss.^(-gamma2-1)./rigidfunc(0.06,U_ss) * (1-0.06);
    
    if do_fullins 
        Cu_ss = 1 - U_ss - G_ss; 
        Ce_ss = 1 - U_ss - G_ss; 
    end
end

%% Set the Grid for Endogenous and Exogenous Variables
LowerBound = [Umin -2*sigma_epsg/sqrt(1-rho_g^2) -2*sigma_epsz/sqrt(1-rho_z^2)];
UpperBound = [Umax 2*sigma_epsg/sqrt(1-rho_g^2) 2*sigma_epsz/sqrt(1-rho_z^2)];

fspace = fundefn(ApproxType,Order,LowerBound,UpperBound,splineorder);
[Grid,nodes]  = funnode(fspace);

% Quadrature Points and Weights to compute expectations wrt supply and demand
[QuadrPoints,QuadrWeights] = qnwnorm([nQuadr nQuadr],[-sigma_epsg^2/2;-sigma_epsz^2/2],[sigma_epsg^2 0; 0 sigma_epsz^2]);

% Build the grid of possible realizations of future shocks
gGrid   = BigGrid(nQuadr^nshocks, Grid(:,2));
zGrid   = BigGrid(nQuadr^nshocks, Grid(:,3));
gNext   = min(max(rho_g*gGrid+kron(ones(size(Grid,1),1),QuadrPoints(:,1)),LowerBound(2)),UpperBound(2));
zNext   = min(max(rho_z*zGrid+kron(ones(size(Grid,1),1),QuadrPoints(:,2)),LowerBound(3)),UpperBound(3));

ExoBigGrid = [gNext zNext];
ExoProbMat = QuadrWeights';

%% Set the initial guess for parameters
% par0 is initial vector of parameters

DepVarX   = log(ones(size(Grid,1),1)*((1-prob_eu).*Ce_ss + prob_eu*theta).^(-sig));
DepVar_Unempf = U_ss*ones(size(Grid,1),1);
DepVar = [DepVarX DepVar_Unempf];

% generate basis functions Chebbasis at Grid where DepVar was defined:
Chebbasis = funbase(fspace,Grid);
% get starting values of parameters by least squares:
par0 = Chebbasis\DepVar;
par0 = par0(:);

%% SOLVE the model
option=optimset('MaxFunEvals',25000,'Display','iter', 'TolFun',1e-12);
par  = fsolve(@model_colloc,par0,option,Grid,fspace,ExoBigGrid,ExoProbMat);

%% SIMULATE the model
t_simul = 500;                            
n_simul = 10;
epsz_Simul = -sigma_epsz^2/2 + sigma_epsz*randn(t_simul,n_simul);
epsg_Simul = -sigma_epsg^2/2 + sigma_epsg*randn(t_simul,n_simul);

for j=1:n_simul
    g_Simul_lag = 0;
    z_Simul_lag = 0;
    Unemp_lag  = U_ss;
    for t=1:t_simul
        z_Simul(t,j) = rho_z*z_Simul_lag + epsz_Simul(t,j);
        g_Simul(t,j) = rho_g*g_Simul_lag + epsg_Simul(t,j);
        Grid_simul = [Unemp_lag g_Simul(t,j) z_Simul(t,j)];
        [Y_simul(t,j),pi_simul(t,j),Unemp_simul(t,j),Ce_simul(t,j)] = GetVar(par,fspace,Grid_simul);
        Unemp_lag = Unemp_simul(t,j);
        
        g_Simul_lag = g_Simul(t,j);
        z_Simul_lag = z_Simul(t,j);
    end
end
G_simul = G_ss*exp(g_Simul);
C_simul = Y_simul - G_ss*exp(g_Simul);
U_mean = mean(mean(Unemp_simul));

%% Calculate multipliers in booms and recessions
t_IRF=20;
g_IRF = [sigma_epsg*rho_g.^(0:t_IRF)'];
griditer = [-2:0.1:0];

for j=1:length(griditer)
    U_ini = U_ss;
    z_IRF(:,j) = griditer(j)*sigma_epsz*rho_z.^(0:t_IRF)';
    
    %U_ini = U_ss+griditer(j);
    %z_IRF(:,j) = zeros(t_IRF+1,1);
    
    Unemp_lag_exp = U_ini; Unemp_lag_contr = U_ini; Unemp_lag_base = U_ini;
    
    for t=1:t_IRF+1
        Grid_IRF_exp = [Unemp_lag_exp g_IRF(t) z_IRF(t,j)];
        Grid_IRF_contr = [Unemp_lag_contr -g_IRF(t) z_IRF(t,j)];
        Grid_IRF_base = [Unemp_lag_base 0 z_IRF(t,j)];
        
        [Y_IRF_base(t,j),pi_IRF_base(t,j),Unemp_IRF_base(t,j),Ce_IRF_base(t,j)] = GetVar(par,fspace,Grid_IRF_base);
        [Y_IRF_exp(t,j),pi_IRF_exp(t,j),Unemp_IRF_exp(t,j), Ce_IRF_exp(t,j)] = GetVar(par,fspace,Grid_IRF_exp);
        [Y_IRF_contr(t,j),pi_IRF_contr(t,j),Unemp_IRF_contr(t,j), Ce_IRF_contr(t,j)] = GetVar(par,fspace,Grid_IRF_contr);
        Unemp_lag_exp = Unemp_IRF_exp(t,j);
        Unemp_lag_contr = Unemp_IRF_contr(t,j);
        Unemp_lag_base = Unemp_IRF_base(t,j);
    end
    
    mult_exp(j) = sum(Y_IRF_exp(:,j)-Y_IRF_base(:,j))/sum(G_ss*(exp(g_IRF)-1));
    mult_contr(j) = sum(Y_IRF_contr(:,j)-Y_IRF_base(:,j))/sum(G_ss*(exp(-g_IRF)-1));
    Ceshare(j)=mean(Ce_IRF_base(:,j)./(1-Unemp_IRF_base(:,j)-G_ss));    
    
end
%% store results
U_ref = U_mean; 
if do_constantelast || do_fullins; U_ref=0.06; end

INDICATOR_Model = Unemp_IRF_base(1,:) - U_ref;


multipliers_for_table=[ mult_exp(:) mult_contr(:)];
indicators=[6 18];
temp=flipud(multipliers_for_table(indicators,:));
table4_entries=temp(:);
results_main_dsge_ce=table4_entries';

save results_main_ce results_main_dsge_ce


clc
close all;
disp('constant elasticity:')
table4_entries'
if ~do_constantelast && ~do_fullins
save sumstats_Model.mat INDICATOR_Model mult_exp mult_contr
do_FigureMultiplier;
end

