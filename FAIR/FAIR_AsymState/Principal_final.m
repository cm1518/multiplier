% Main prg to run  estimation with asymmetry and state dependence.

tic

clear all
close all
clc


%seed

load seed_for_run 

rng(s)

%setup structure that encodes prior etc

load setup_for_state_asym


%% estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );

save state_asym_results
toc

Figure_6;
