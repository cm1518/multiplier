% Main prg to run the symmetric estimation

clear all
close all
clc




load seed_for_run 

rng(s)

load setup_sym


%estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );



save sym_paper



compute_multiplier
