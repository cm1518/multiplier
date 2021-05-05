% Main prg to run the symmetric estimation

clear all
close all
clc
tic

%seed

load seed_for_run
rng(s);


load setup_sym



[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );
save SUR_sym_paper
toc


multiplier;