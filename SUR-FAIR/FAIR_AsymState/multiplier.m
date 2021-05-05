
lower=5;
upper=95;

setup.number_of_draws=setup.keep_draw*size(draws,2);
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,1000,1);

YGr=1; 
 
hor=20;

%indicator:
ind=0;

IRF_neg=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));
IRF_pos=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));


Msum=[];Mmax=[];
Ndraws=1000;

index=0;

 
INDICATOR=[-1 2];

for ind=INDICATOR
    index=index+1;
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
            
            
            gamma_neg{jj}=params(counter);
            counter=counter+1;
            
            beta_neg{jj}=params(counter);
            counter=counter+1;
            
            gamma_pos{jj}=params(counter);
            counter=counter+1;
            beta_pos{jj}=params(counter);
            counter=counter+1;
            
        end
        
        
        for jj=1:setup.dep_variables
            Sigma_lagged_neg{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
            for gg=1:setup.num_Gaussian{jj}
                [ Sigma_lagged_neg{jj}(gg,:) ] = NL_IRF(NL_a_neg{jj}(gg),NL_b_neg{jj}(gg),NL_c_neg{jj}(gg),gamma_neg{jj},beta_neg{jj},setup.lags_shocks,ind);
                
            end
            
            Sigma_lagged_pos{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
            for gg=1:setup.num_Gaussian{jj}
                [ Sigma_lagged_pos{jj}(gg,:) ] = NL_IRF(NL_a_pos{jj}(gg),NL_b_pos{jj}(gg),NL_c_pos{jj}(gg),gamma_pos{jj},beta_pos{jj},setup.lags_shocks,ind);
                
            end
            
            
            IRF_neg(jj,:,kk)=sum(Sigma_lagged_neg{jj},1);
            IRF_pos(jj,:,kk)=sum(Sigma_lagged_pos{jj},1);
        end
        
        % Store multipliers matrix, first column with neg shocks, second column for pos shocks
        [mx_pos i_mxpos]=max(YGr*abs(IRF_pos(2,:,kk)));
        [mx_neg i_mxneg]=max(YGr*abs(IRF_neg(2,:,kk)));
        
        Mmax(kk,:,index)=[sign(IRF_neg(2,i_mxneg,kk)).*mx_neg./max(IRF_neg(1,:,kk)) sign(IRF_pos(2,i_mxpos,kk)).*mx_pos./max(IRF_pos(1,:,kk))];
        Msum(kk,:,index)=[sum(YGr*IRF_neg(2,1:hor,kk))./sum(IRF_neg(1,1:hor,kk)) sum(YGr*IRF_pos(2,1:hor,kk))./sum(IRF_pos(1,1:hor,kk))];
        
        
    end
    
    
end

msum=squeeze(median(Msum,1))';
msumlow=squeeze(prctile(Msum,lower,1))';
msumhigh=squeeze(prctile(Msum,upper,1))';

mmax=squeeze(median(Mmax,1))';
mmaxlow=squeeze(prctile(Mmax,lower,1))';
mmaxhigh=squeeze(prctile(Mmax,upper,1))';



%"sum" multiplier, comparing with linear GMA estimate
display('  ')
display('  ')
horM=20;
display(['Hor=',num2str(horM),'       Neg         Pos'])
display(['UR       multiplier msum  '])
display(num2str([INDICATOR' squeeze(prctile(Msum,50))']))

display(['5% lower band:   ']);
display(num2str([INDICATOR' squeeze(prctile(Msum,5))']))

display(['95% upper band:  ']);
display(num2str([INDICATOR' squeeze(prctile(Msum,95))']))


results_asym_state_sur.med=prctile(Msum,50);
results_asym_state_sur.lower=prctile(Msum,5);
results_asym_state_sur.upper=prctile(Msum,95);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m(U high)>m(U low)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Msum

hor_Msum=20;

% Msum(draw#,sign,indicator) 

%Msum_pos

% Count # of points above 45degree line:
count=sum(squeeze(Msum(:,2,2))>squeeze(Msum(:,2,1)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(m_pos(U high)>m_pos(U low) at hor=',num2str(hor_Msum),') = ',num2str(count/length(Msum))])
disp(' ')


results_asym_state_sur.prob_pos=num2str(count/length(Msum));

% Count # of points above 45degree line:
count=sum(squeeze(Msum(:,1,2))>squeeze(Msum(:,1,1)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(m_neg(U high)>m_neg(U low) at hor=',num2str(hor_Msum),') = ',num2str(count/length(Msum))])
disp(' ')

results_asym_state_sur.prob_neg=num2str(count/length(Msum));

save asym_state_results_for_table_sur results_asym_state_sur


