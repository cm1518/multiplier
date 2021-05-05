

setup.number_of_draws=setup.keep_draw*length(draws);

lower=5;
upper=95;

Msum=[];Mmax=[];
YGr=1; 
Ndraws=400;

hor=20;

index1=0;
std_ind=std(setup.indicator);
mu_ind=mean(setup.indicator);
INDICATOR=-1:.5:2;



    
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
            % ind_vec=[zeros(1,setup.lags)];
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

msum=squeeze(median(Msum,3));
msumlow=squeeze(prctile(Msum,lower,3));
msumhigh=squeeze(prctile(Msum,upper,3));

mmax=squeeze(median(Mmax,3));
mmaxlow=squeeze(prctile(Mmax,lower,3));
mmaxhigh=squeeze(prctile(Mmax,upper,3));

disp('median multiplier, negative shock top row, positive shock bottom row, U low left column, U high right column')
msum(:,[1 end])
disp('5th percentile multiplier, negative shock top row, positive shock bottom row, U low left column, U high right column')
msumlow(:,[1 end])
disp('95th percentile multiplier, negative shock top row, positive shock bottom row, U low left column, U high right column')
msumhigh(:,[1 end])

%VAR-based multipliers:
for jj=1:setup.size_obs
    irf_var(:,jj)=squeeze(setup.store_responses(jj,setup.index_unrestricted,:));
end
msum_var=sum(irf_var(1:21,4))./sum(irf_var(1:21,2))*YGr;
mmax_var=max(irf_var(1:21,4))./max(irf_var(1:21,2))*YGr;

EZ=0;



save('SumStats.mat','Msum','Mmax','INDICATOR')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m(U high)>m(U low)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Msum

hor_Msum=20;

% Msum(sign,indicator,draw#) 

%Msum_pos

% Count # of points above 45degree line:
count=sum(squeeze(Msum(2,end,:))>squeeze(Msum(2,1,:)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(m_pos(U high)>m_pos(U low) at hor=',num2str(hor_Msum),') = ',num2str(round(count/length(Msum),1))])
disp(' ')


%Msum_neg

% Count # of points above 45degree line:
count=sum(squeeze(Msum(1,end,:))>squeeze(Msum(1,1,:)));

% Posterior proba that msum(UR high)>msum(UR low)
disp(' ')
disp(['P(m_neg(U high)>m_neg(U low) at hor=',num2str(hor_Msum),') = ',num2str(round(count/length(Msum),1))])
disp(' ')


