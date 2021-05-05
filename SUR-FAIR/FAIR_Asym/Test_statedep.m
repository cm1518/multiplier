%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m_neg>m_pos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters of plot:
ctr_N=10; %density of contour plot

%%%%%%%%%%
hor_Msum=20;


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(:,2,hor_Msum))>squeeze(MsumS(:,1,hor_Msum)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(msum_neg>msum_pos at hor=',num2str(hor_Msum),') = ',num2str(1-count/length(MsumS))])
disp(' ')

results_asym_sur.prob=1-count/length(MsumS);

save asym_results_for_table_sur results_asym_sur
