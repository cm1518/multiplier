% Main prg to run the BM asymmetric estimation

clear all
close all
clc
tic

%seed
load seed_for_run
rng(s);



load setup_asym

[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );
save SUR_asym_paper
toc

multiplier;
