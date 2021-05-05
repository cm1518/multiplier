% estimates some of the values for table 3

tic

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load setup_for_asym_state_47;

%seed



load seed_for_run;
rng(s);

%% estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );

save asym_state_47_results;
toc



Test_multiplier_statedep;
