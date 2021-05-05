%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m_neg>m_pos for Msum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters of plot:
ctr_N=10; %density of contour plot


%%%%%
hor_Msum=20;


% Count # of points above 45degree line:
count=sum(squeeze(MsumS(1,:,hor_Msum))>squeeze(MsumS(2,:,hor_Msum)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(msum_neg>msum_pos at hor=',num2str(hor_Msum),') = ',num2str(1-count/length(MsumS))])
disp(' ')








% Value of multiplier and error-bands

%horM=20;
horM=20;
display('  ')
display('  ')
display(['Hor=',num2str(horM),'       Pos         Neg'])
temp=prctile(MsumS,50,2);
mult_med=temp(:,:,horM);
display(['median  ',num2str(squeeze(temp(:,:,horM))')]);

temp=prctile(MsumS,5,2);
mult_lower=temp(:,:,horM);
display(['5%  ',num2str(squeeze(temp(:,:,horM))')]);
temp=prctile(MsumS,95,2);
mult_upper=temp(:,:,horM);
display(['95%   ',num2str(squeeze(temp(:,:,horM))')]);
temp=std(MsumS,0,2);
display(['std   ',num2str(squeeze(temp(:,:,horM))')]);


results_asym.lower=mult_lower;
results_asym.med=mult_med;
results_asym.upper=mult_upper;
results_asym.prob=1-count/length(MsumS);

save asym_results_for_table results_asym
