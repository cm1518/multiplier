% Main prg to run the BM asymmetric estimation

clear all
close all
clc
tic

load seed_for_run
rng(s)



load setup_asym_state


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );

save SUR_asym_state_paper

multiplier;



toc
