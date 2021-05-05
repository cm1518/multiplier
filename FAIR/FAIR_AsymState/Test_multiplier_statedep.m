
setup.number_of_draws=setup.keep_draw*length(draws);

lower=5;
upper=95;

hor=20;

Msum=[];Mmax=[];
YGr=1;
Ndraws=100;

index1=0;



load seed_for_table
rng(s)


INDICATOR=[-1 2];



    
    IRFs=[];
    for sgn=-1:2:1
        index1=index1+1;
        
        index2=0;
        for indicator=INDICATOR
            index2=index2+1;
            inds=setup.number_of_draws/setup.keep_draw;
            indices_for_draws=unidrnd(inds,Ndraws,1);
            
            mean_ind=0; %adjust if indicator mean is not 0
            for kk=1:length(indices_for_draws)
                size_of_shock=sgn;
               
                shock_contemp=zeros(setup.size_obs,1);
                shock_contemp(setup.index_unrestricted)=size_of_shock;
                
                ind_vec=indicator*[ones(1,setup.lags+1)];
                for jj=1:setup.lags+1
                    
                    
                    epsilon=[zeros(setup.size_obs,jj-1) shock_contemp zeros(setup.size_obs,setup.lags+1-jj)];
                    [ Sigma, intercept] = unwrap_NL_IRF( draws(:,indices_for_draws(kk)),epsilon,setup ,ind_vec,mean_ind);
                  
                    
                    IRFs(:,jj,kk)=Sigma(:,setup.index_unrestricted,jj);
                end
                
                 
                Msum(index1,index2,kk)=sum(YGr*IRFs(4,1:hor,kk),2)./sum(IRFs(2,1:hor,kk),2);
                Mmax(index1,index2,kk)=max(YGr*IRFs(4,1:hor,kk))./max(IRFs(2,1:hor,kk));
                
                clear ind_vec
                
            end
            
            IRFsmed=squeeze(median(IRFs,3));
            IRFslow=squeeze(prctile(IRFs,lower,3));
            IRFshigh=squeeze(prctile(IRFs,upper,3));
            
               
            
        end
    end
    save Test_statedep.mat Msum INDICATOR
    


msum=squeeze(median(Msum,3));
msumlow=squeeze(prctile(Msum,lower,3));
msumhigh=squeeze(prctile(Msum,upper,3));

mmax=squeeze(median(Mmax,3));
mmaxlow=squeeze(prctile(Mmax,lower,3));
mmaxhigh=squeeze(prctile(Mmax,upper,3));


for jj=1:setup.size_obs
    irf_var(:,jj)=squeeze(setup.store_responses(jj,setup.index_unrestricted,:));
end
msum_var=sum(irf_var(1:21,4))./sum(irf_var(1:21,2))*YGr;
mmax_var=max(irf_var(1:21,4))./max(irf_var(1:21,2))*YGr;

EZ=0;




display('  ')
display('  ')
horM=20;
display(['Hor=',num2str(horM),'       Neg         Pos'])
display(['UR       multiplier msum  '])
display(num2str([INDICATOR' prctile(Msum,50,3)']))

display(['5% lower band:   ']);
display(num2str([INDICATOR' prctile(Msum,5,3)']))

display(['95% upper band:  ']);
display(num2str([INDICATOR' prctile(Msum,95,3)']))

results_asym_state.med=prctile(Msum,50,3);
results_asym_state.lower=prctile(Msum,5,3);
results_asym_state.upper=prctile(Msum,95,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m(U high)>m(U low)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Msum

hor_Msum=20;




% Count # of points above 45degree line:
count=sum(squeeze(Msum(2,2,:))>squeeze(Msum(2,1,:)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(m_pos(U high)>m_pos(U low) at hor=',num2str(hor_Msum),') = ',num2str(count/length(Msum))])
disp(' ')

results_asym_state.prob_pos=num2str(count/length(Msum));
% Count # of points above 45degree line:
count=sum(squeeze(Msum(1,2,:))>squeeze(Msum(1,1,:)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(m_neg(U high)>m_neg(U low) at hor=',num2str(hor_Msum),') = ',num2str(count/length(Msum))])
disp(' ')
results_asym_state.prob_neg=num2str(count/length(Msum));

save asym_state_results_for_table results_asym_state
