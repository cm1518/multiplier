% Main prg to run the BDM asymmetric estimation

clear all
close all
clc

%set seed
rng(1)

tic


%load setup structure

load setup_for_asym;

%estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );





toc


save results_asym


Figure_2


