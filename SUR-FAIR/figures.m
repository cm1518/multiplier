
%Figures involving SUR approach

clear all
close all
clc


load seed_for_figures 
rng(s);


addpath('..')

set(0,'defaulttextInterpreter','latex')

cd FAIR_Sym
plots_irfs_paper;

cd ../FAIR_Asym
addpath('..')
clear all

load SUR_asym_paper

load ('../FAIR_Sym/SymResults.mat')
plots_irfs_paper

cd ..
%Get Asymmetric IRFs from state dep code (and load previous results to get graph for paper)
cd FAIR_AsymState
clear all
load SUR_asym_state_paper
load ('../FAIR_Sym/SymResults.mat')

Plots_st
PlotIRFforPaper_shadedbands_Revision

cd ..





