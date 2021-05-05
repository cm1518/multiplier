%computes multipliers and stores estimated IRFs for later plotting


setup.number_of_draws=setup.keep_draw*size(draws,2);

inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,5000,1);

YGr=1;
hor=30;
indicator=0;

IRF_neg=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));
for kk=1:length(indices_for_draws)
    params=draws(:,indices_for_draws(kk));
    counter=1;
    for jj=1:setup.dep_variables
        
        %unwrapping parameters
        
        intercept{jj}=params(counter);
        counter=1+counter;
        
        
        
        lagged_y_param{jj}=params(counter:counter+setup.lags_endo-1);
        counter=counter+setup.lags_endo;
        
        lagged_exo_param{jj}=params(counter:counter+setup.lags_exo*setup.dim_exo-1); %I assume that parameters here are ordered by variable first, then by lag
        counter=counter+setup.lags_exo*setup.dim_exo;
        
        
        if jj==1
            counter_linear=counter; %number of 'linear' (except AR coefficients) plus 1, same for all dep. variables
        end
        
        NL_a_neg{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_b_neg{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_c_neg{jj}=zeros(setup.num_Gaussian{jj},1);
        
        
        for ll=1:setup.num_Gaussian{jj}
            NL_a_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_b_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_c_neg{jj}(ll)=params(counter);
            counter=counter+1;
        end
        
        
    end
    
    
    for jj=1:setup.dep_variables
        Sigma_lagged_neg{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
        for gg=1:setup.num_Gaussian{jj}
            [ Sigma_lagged_neg{jj}(gg,:) ] = NL_IRF(NL_a_neg{jj}(gg),NL_b_neg{jj}(gg),NL_c_neg{jj}(gg),setup.lags_shocks);
            
        end
        
        IRF_neg(jj,:,kk)=(sum(Sigma_lagged_neg{jj},1));
    end
    
    Mmax(kk,:)=[max(YGr*IRF_neg(2,:,kk))./max(IRF_neg(1,:,kk))];
    Msum(kk,:)=[sum(YGr*IRF_neg(2,1:hor,kk))./sum(IRF_neg(1,1:hor,kk))];
    
    for h=1:hor
        MsumS(kk,h)=[sum(YGr*IRF_neg(2,1:h,kk))./sum(IRF_neg(1,1:h,kk)) ];
    end
    
end

%moments of irfs
upper=95;
lower=5;

median_neg=prctile(IRF_neg,50,3); %first dimension is the observable, second the horizon
lower_neg=prctile(IRF_neg,lower,3);
upper_neg=prctile(IRF_neg,upper,3);


mmax=squeeze(median(Mmax,1));
mmaxlow=squeeze(prctile(Mmax,lower,1));
mmaxhigh=squeeze(prctile(Mmax,upper,1));

msum=squeeze(median(Msum,1));
msumlow=squeeze(prctile(Msum,lower,1));
msumhigh=squeeze(prctile(Msum,upper,1));





%moments of irfs
upper=95;
lower=5;

mmax=[squeeze(prctile(Mmax,50,1))];
mmaxlow=[squeeze(prctile(Mmax,lower,1))];
mmaxhigh=[squeeze(prctile(Mmax,upper,1))];

msum=squeeze(prctile(Msum,50,1));
msumlow=squeeze(prctile(Msum,lower,1));
msumhigh=squeeze(prctile(Msum,upper,1));


%m_sum:
msum=prctile(MsumS,50,1);
msum_ub=(prctile(MsumS,upper,1));
msum_lb=(prctile(MsumS,lower,1));

MsumHor=[];
MsumHor=squeeze(msum)';
MsumHor_lb=squeeze(msum_lb)';
MsumHor_ub=squeeze(msum_ub)';







%Store symmetric results
irfsym=prctile(IRF_neg,50,3); %first dimension is the observable, second the horizon
irfsym_lower=prctile(IRF_neg,lower,3);
irfsym_upper=prctile(IRF_neg,upper,3);



display('  ')
display('  ')
horM=20;
display(['Hor=',num2str(horM),'             '])
display(['multiplier  ', num2str([ prctile(MsumS(:,horM),50)])])
display(['5% lower band:   ', num2str([ prctile(MsumS(:,horM),5)])]);
display(['95% upper band:  ', num2str([ prctile(MsumS(:,horM),95)])]);

save('SymResults.mat','irfsym','irfsym_upper','irfsym_lower')

results_sur_sym.lower=prctile(MsumS(:,horM),5);
results_sur_sym.med=prctile(MsumS(:,horM),50);
results_sur_sym.upper=prctile(MsumS(:,horM),95);


save sym_results_for_table_sur results_sur_sym



cd ..
save('SymResults.mat','irfsym','irfsym_upper','irfsym_lower')



