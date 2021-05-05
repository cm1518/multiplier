% Main prg to run the BM asymmetric estimation

clear all
clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Part to change depending on whether RB or CM is running the code
fileroot='E:\BM_non_linear\single_equation_asymmetry';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization/Parametrization of routine
setup_single_equation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find initial values via OLS
data=load(setup.data);
data=data.data;
setup.sample_size=size(data,2);
setup.dim_add_matrices=[setup.sample_size 1];
OLS_estimate=((data(2:end,:))*data(2:end,:)')\(data(2:end,:))*data(1,:)';

%match OLS estimates of shock responses to a, b and c coefficients
options = optimset('Display','iter','TolFun',1e-20,'TolX',1e-20,'MaxIter',5000, 'MaxFunEvals',500000);
opt=@(params)IRF_match( params,setup,OLS_estimate(2+setup.lags_endo+setup.lags_exo*setup.dim_exo:end)) ;
[xestimate,functionvalue1]=fminsearch(opt,ones(3,1),options);

setup.initial_parameter=[OLS_estimate(1:1+setup.lags_endo+setup.lags_exo*setup.dim_exo);xestimate;xestimate;1];
setup.length_param_vector=length(setup.initial_parameter);
%estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );

%add_matrices are the estimated shocks






save test_single_eq
toc
