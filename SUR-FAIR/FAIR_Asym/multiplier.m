

setup.number_of_draws=setup.keep_draw*size(draws,2);
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,5000,1);

YGr=1;

hor=30;

indicator=0;



IRF_neg=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));
IRF_pos=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));

Mmax=[]; Msum=[];

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
        
        
        
        NL_a_pos{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_b_pos{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_c_pos{jj}=zeros(setup.num_Gaussian{jj},1);
        
        for ll=1:setup.num_Gaussian{jj}
            NL_a_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_b_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_c_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_a_pos{jj}(ll)=params(counter);
            counter=counter+1;
            NL_b_pos{jj}(ll)=params(counter);
            counter=counter+1;
            NL_c_pos{jj}(ll)=params(counter);
            counter=counter+1;
        end
        
        
    end
    
    
    for jj=1:setup.dep_variables
        Sigma_lagged_neg{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
        for gg=1:setup.num_Gaussian{jj}
            [ Sigma_lagged_neg{jj}(gg,:) ] = NL_IRF(NL_a_neg{jj}(gg),NL_b_neg{jj}(gg),NL_c_neg{jj}(gg),setup.lags_shocks);
            
        end
        
        Sigma_lagged_pos{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
        for gg=1:setup.num_Gaussian{jj}
            [ Sigma_lagged_pos{jj}(gg,:) ] = NL_IRF(NL_a_pos{jj}(gg),NL_b_pos{jj}(gg),NL_c_pos{jj}(gg),setup.lags_shocks);
            
        end
        
        
        IRF_neg(jj,:,kk)=(sum(Sigma_lagged_neg{jj},1));
        IRF_pos(jj,:,kk)=sum(Sigma_lagged_pos{jj},1);
    end
    
    % Store multipliers matrix, first column with neg shocks, second column for pos shocks
    [mx_pos i_mxpos]=max(YGr*abs(IRF_pos(2,:,kk)));
    [mx_neg i_mxneg]=max(YGr*abs(IRF_neg(2,:,kk)));
    
    Mmax(kk,:)=[sign(IRF_neg(2,i_mxneg,kk)).*mx_neg./max(IRF_neg(1,:,kk)) sign(IRF_pos(2,i_mxpos,kk)).*mx_pos./max(IRF_pos(1,:,kk))];
    Msum(kk,:)=[sum(YGr*IRF_neg(2,:,kk))./sum(IRF_neg(1,:,kk)) sum(YGr*IRF_pos(2,:,kk))./sum(IRF_pos(1,:,kk))];
    
    for h=1:hor
        MsumS(kk,:,h)=[sum(YGr*IRF_neg(2,1:h,kk))./sum(IRF_neg(1,1:h,kk)) sum(YGr*IRF_pos(2,1:h,kk))./sum(IRF_pos(1,1:h,kk))];
    end
end

%moments of irfs
upper=97.5;
lower=2.5;

median_neg=prctile(IRF_neg,50,3); %first dimension is the observable, second the horizon
lower_neg=prctile(IRF_neg,lower,3);
upper_neg=prctile(IRF_neg,upper,3);

median_pos=prctile(IRF_pos,50,3);
lower_pos=prctile(IRF_pos,lower,3);
upper_pos=prctile(IRF_pos,upper,3);







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

MsumHor=squeeze(msum)';
MsumHor_lb=squeeze(msum_lb)';
MsumHor_ub=squeeze(msum_ub)';




%"sum" multiplier, comparing with linear FAIR estimate
display('  ')
display('  ')
horM=20;
display(['Hor=',num2str(horM),'             Neg          Pos'])
display(['multiplier msum  ', num2str([ prctile(MsumS(:,:,horM),50)])])
display(['5% lower band:   ', num2str([ prctile(MsumS(:,:,horM),5)])]);
display(['95% upper band:  ', num2str([ prctile(MsumS(:,:,horM),95)])]);
display(['Std dev:         ', num2str([ std(MsumS(:,:,horM)) ])]);

results_asym_sur.lower=prctile(MsumS(:,:,horM),5);
results_asym_sur.med=prctile(MsumS(:,:,horM),50);
results_asym_sur.upper=prctile(MsumS(:,:,horM),95);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m_neg>m_pos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Test_statedep
