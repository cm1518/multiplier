clear;
clc;
close all

%collect results for table 1

load FAIR_Sym/sym_results_for_table.mat

load FAIR_Asym/asym_results_for_table.mat

load ../SUR-FAIR/FAIR_Sym/sym_results_for_table_sur.mat

load ../SUR-FAIR/FAIR_Asym/asym_results_for_table_sur.mat

lower=round([results_sym.lower*ones(2,1)';results_asym.lower';results_sur_sym.lower*ones(2,1)';fliplr(results_asym_sur.lower)],1);
med=round([results_sym.med*ones(2,1)';results_asym.med';results_sur_sym.med*ones(2,1)';fliplr(results_asym_sur.med)],2);
upper=round([results_sym.upper*ones(2,1)';results_asym.upper';results_sur_sym.upper*ones(2,1)';fliplr(results_asym_sur.upper)],1);
prob=round([0;results_asym.prob;0;results_asym_sur.prob],2);
specification={'symmetric rec. FAIR';'asymmetric rec. FAIR';'symmetric narr. FAIR';'asymmetric narr. FAIR'};
lower_expansionary=lower(:,1);
lower_contractionary=lower(:,2);
upper_expansionary=upper(:,1);
upper_contractionary=upper(:,2);
med_expansionary=med(:,1);
med_contractionary=med(:,2);


T = table(specification,lower_expansionary,lower_contractionary,med_expansionary,med_contractionary,upper_expansionary,upper_contractionary,prob);
writetable(T,'table_1.xlsx');

%table 2

load FAIR_AsymState/asym_state_results_for_table.mat
load ../SUR-FAIR/FAIR_AsymState/asym_state_results_for_table_sur.mat

lower=round([results_asym_state.lower;squeeze(results_asym_state_sur.lower)],1);
upper=round([results_asym_state.upper;squeeze(results_asym_state_sur.upper)],1);
med=round([results_asym_state.med;squeeze(results_asym_state_sur.med)],2);

rec_med=med([2 6 1 5]);
narr_med=med([4 8 3 7]);

rec_lower=lower([2 6 1 5]);
narr_lower=lower([4 8 3 7]);

rec_upper=upper([2 6 1 5]);
narr_upper=upper([4 8 3 7]);

probs_rec=round([str2num(results_asym_state.prob_pos) str2num(results_asym_state.prob_pos) str2num(results_asym_state.prob_neg) str2num(results_asym_state.prob_neg)],3);
probs_narr=round([str2num(results_asym_state_sur.prob_pos) str2num(results_asym_state_sur.prob_pos) str2num(results_asym_state_sur.prob_neg) str2num(results_asym_state_sur.prob_neg)],3);


low_expansionary=[rec_med(1);narr_med(1);rec_lower(1);narr_lower(1);rec_upper(1);narr_upper(1);probs_rec(1);probs_narr(1)];
high_expansionary=[rec_med(2);narr_med(2);rec_lower(2);narr_lower(2);rec_upper(2);narr_upper(2);probs_rec(2);probs_narr(2)];

low_contractionary=[rec_med(3);narr_med(3);rec_lower(3);narr_lower(3);rec_upper(3);narr_upper(3);probs_rec(3);probs_narr(3)];
high_contractionary=[rec_med(4);narr_med(4);rec_lower(4);narr_lower(4);rec_upper(4);narr_upper(4);probs_rec(4);probs_narr(4)];

specification={'asymmetric rec. FAIR with state dep. median';'asymmetric narr. FAIR with state dep. median';'asymmetric rec. FAIR with state dep. 5th perc.';'asymmetric narr. FAIR with state dep. 5th perc.';'asymmetric rec. FAIR with state dep. 95th perc.';'asymmetric narr. FAIR with state dep. 95th perc.';'probability that multiplier with U high is larger than with U low, rec.';'probability that multiplier with U high is larger than with U low, narr.'};
T2=table(specification,low_expansionary,high_expansionary,low_contractionary,high_contractionary);

writetable(T2,'table_2.xlsx');

%table 3


load FAIR_AsymState_47/asym_state_results_for_table_47.mat
load ../SUR-FAIR/FAIR_AsymState_47/asym_state_results_for_table_sur_47.mat

lower=round([results_asym_state_47.lower;squeeze(results_asym_state_sur_47.lower)],1);
upper=round([results_asym_state_47.upper;squeeze(results_asym_state_sur_47.upper)],1);
med=round([results_asym_state_47.med;squeeze(results_asym_state_sur_47.med)],2);

rec_med=med([2 6 1 5]);
narr_med=med([4 8 3 7]);

rec_lower=lower([2 6 1 5]);
narr_lower=lower([4 8 3 7]);

rec_upper=upper([2 6 1 5]);
narr_upper=upper([4 8 3 7]);

probs_rec=round([str2num(results_asym_state_47.prob_pos) str2num(results_asym_state_47.prob_pos) str2num(results_asym_state_47.prob_neg) str2num(results_asym_state_47.prob_neg)],3);
probs_narr=round([str2num(results_asym_state_sur_47.prob_pos) str2num(results_asym_state_sur_47.prob_pos) str2num(results_asym_state_sur_47.prob_neg) str2num(results_asym_state_sur_47.prob_neg)],3);


low_expansionary=[rec_med(1);narr_med(1);rec_lower(1);narr_lower(1);rec_upper(1);narr_upper(1);probs_rec(1);probs_narr(1)];
high_expansionary=[rec_med(2);narr_med(2);rec_lower(2);narr_lower(2);rec_upper(2);narr_upper(2);probs_rec(2);probs_narr(2)];

low_contractionary=[rec_med(3);narr_med(3);rec_lower(3);narr_lower(3);rec_upper(3);narr_upper(3);probs_rec(3);probs_narr(3)];
high_contractionary=[rec_med(4);narr_med(4);rec_lower(4);narr_lower(4);rec_upper(4);narr_upper(4);probs_rec(4);probs_narr(4)];

specification={'asymmetric rec. FAIR with state dep. median';'asymmetric narr. FAIR with state dep. median';'asymmetric rec. FAIR with state dep. 5th perc.';'asymmetric narr. FAIR with state dep. 5th perc.';'asymmetric rec. FAIR with state dep. 95th perc.';'asymmetric narr. FAIR with state dep. 95th perc.';'probability that multiplier with U high is larger than with U low, rec.';'probability that multiplier with U high is larger than with U low, narr.'};
T3=table(specification,low_expansionary,high_expansionary,low_contractionary,high_contractionary);

writetable(T3,'table_3.xlsx');

%table 4        

load ../theoretical_model/results_main
load ../theoretical_model/results_main_pi
load ../theoretical_model/results_main_ce

specification={'Baseline'; 'Perfect Insurance';'Constant Elasticity'};

expansionary_low=round([results_main_dsge(1); results_main_dsge_pi(1);results_main_dsge_ce(1)],2);

expansionary_high=round([results_main_dsge(2); results_main_dsge_pi(2);results_main_dsge_ce(2)],2);

contractionary_low=round([results_main_dsge(3); results_main_dsge_pi(3);results_main_dsge_ce(3)],2);

contractionary_high=round([results_main_dsge(4); results_main_dsge_pi(4);results_main_dsge_ce(4)],2);

T4=table(specification,expansionary_low,expansionary_high,contractionary_low,contractionary_high);

writetable(T4,'table_4.xlsx');
