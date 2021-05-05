function [Y,pi,Unemp,Ce,x,resid_f,G,a,z, Cu] = GetVar(par,fspace,Grid)
global phipi sig Cu_ss Ce_ss U_ss G_ss do_fullins do_constantelast

FUNCS     = ComputePolFunc(par,fspace,Grid); %expectations of RHS of Euler equation
x = FUNCS(:,1);
Unemp_f = FUNCS(:,2);

FUNCS_ss = ComputePolFunc(par,fspace,[U_ss 0 0]);
x_ss = FUNCS_ss(:,1);

Unemp_lag = Grid(:,1); % lag unemployemnt (not used)
G         = G_ss*exp(Grid(:,2));
z         = Grid(:,3);

a         = 0; % log of TFP (constant)

% value of unemployment consumption
Cu      = Cu_ss*ones(size(Grid,1),1);
if do_fullins; Cu = 1 - U_ss - G; end

% variables if constraint is not binding
Unemp_n       = U_ss;
%Unemp_n     =prob_eu*(1-Unemp_lag) + (1-prob_ue)*Unemp_lag; %natural rate of unemployment (without nominal rigidity)
Ce_n         = 1-(Unemp_n.*Cu+G)./(1-Unemp_n); %consumption of employed worker at the natural rate
rbar         = -sig*log(Ce_ss)-x_ss ; % interest rate such that inflation equal target in steady state
pi_n         = 1/phipi*(-sig*log(Ce_n) + z - x - rbar);
rigid_n      = rigidfunc(Unemp_n,Unemp_n);

isbinding = logical(pi_n < (log(rigid_n) -a));
if do_constantelast; isbinding = true(size(Grid,1),1); end

% variables if constraint is binding
rigid_f =  rigidfunc(Unemp_f,Unemp_n);
pi_f  = (log(rigid_f) -a);
Ce_f  = exp(-1/sig*(phipi*(pi_f)+rbar + x - z));

if do_fullins; Cu(isbinding,1) = Ce_f(isbinding,1); end

resid_f = - Unemp_f  + (1 - Ce_f - G)./(1+Cu - Ce_f); % resource constraint + production function

Ce(isbinding,1) = Ce_f(isbinding);
Ce(~isbinding,1) = Ce_n(~isbinding);
pi(isbinding,1) = pi_f(isbinding);
pi(~isbinding,1) = pi_n(~isbinding);


Unemp  = (1 - Ce - G)./(1+Cu - Ce); 

Y = 1-Unemp;

