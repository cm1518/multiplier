%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot state-dep multiplier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load test_NL_042816

% load test_NL


if exist('SumStats.mat')~=2
    
    multiplier
else
    load SumStats
end


lower=16;
upper=84;

lower2=5;
upper2=95;

msum=squeeze(median(Msum,3));
msumlow=squeeze(prctile(Msum,lower,3));
msumhigh=squeeze(prctile(Msum,upper,3));
msumlow2=squeeze(prctile(Msum,lower2,3));
msumhigh2=squeeze(prctile(Msum,upper2,3));

mmax=squeeze(median(Mmax,3));
mmaxlow=squeeze(prctile(Mmax,lower,3));
mmaxhigh=squeeze(prctile(Mmax,upper,3));
mmaxlow2=squeeze(prctile(Mmax,lower2,3));
mmaxhigh2=squeeze(prctile(Mmax,upper2,3));


plot_shocks
Shock_dist_wrt_indicator

% Plot_statedep_withdist_detrendedUR
Plot_statedep_withdist_UR

% Test state dependence
Test_multiplier_statedep

