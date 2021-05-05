

setup.number_of_draws=setup.keep_draw*size(draws,2);
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,1000,1);

YGr=1; 
%indicator:
ind=-2;

hor=30;

IRF_neg=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));
IRF_pos=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));

Mmax=[]; Msum=[];



index_ind=0;
for indicator=[-1 0 2];
    
    index_ind=index_ind+1;
    
    ind=indicator;
    
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
        
        %Re-scale by size of G shock
        IRFs(2,:,kk,2)=IRF_neg(2,:,kk)/max(abs(IRF_neg(1,:,kk)));
        IRFs(2,:,kk,1)=IRF_pos(2,:,kk)/max(abs(IRF_pos(1,:,kk)));
        
        IRFs(1,:,kk,2)=IRF_neg(1,:,kk);
        IRFs(1,:,kk,1)=IRF_pos(1,:,kk);
                
        % Store multipliers matrix, first column with neg shocks, second column for pos shocks
        [mx_pos i_mxpos]=max(YGr*abs(IRF_pos(2,:,kk)));
        [mx_neg i_mxneg]=max(YGr*abs(IRF_neg(2,:,kk)));
        
        Mmax(kk,:)=[sign(IRF_neg(2,i_mxneg,kk)).*mx_neg./max(IRF_neg(1,:,kk)) sign(IRF_pos(2,i_mxpos,kk)).*mx_pos./max(IRF_pos(1,:,kk))];
        Msum(kk,:)=[sum(YGr*IRF_neg(2,:,kk))./sum(IRF_neg(1,:,kk)) sum(YGr*IRF_pos(2,:,kk))./sum(IRF_pos(1,:,kk))];
        
        
    end
    
  
    IRFsmed=squeeze(prctile(IRFs,50,3));
    IRFslow=squeeze(prctile(IRFs,lower,3));
    IRFshigh=squeeze(prctile(IRFs,upper,3));
    
    
    IRy_med(:,:,index_ind)=squeeze(IRFsmed(2,:,:));
    IRy_lb(:,:,index_ind)=squeeze(IRFslow(2,:,:));
    IRy_ub(:,:,index_ind)=squeeze(IRFshigh(2,:,:));
    
    
    IRg_med(:,:,index_ind)=squeeze(IRFsmed(1,:,:));
    IRg_lb(:,:,index_ind)=squeeze(IRFslow(1,:,:));
    IRg_ub(:,:,index_ind)=squeeze(IRFshigh(1,:,:));
    
end


%Now order the IRFs int he order that they will be plotted:
irf_med(:,1)=IRg_med(:,1,2)/max(IRg_med(:,1,2));
irf_med(:,2:4)=IRy_med(:,1,1:3);
irf_med(:,5)=IRg_med(:,2,2)/max(IRg_med(:,2,2));
irf_med(:,6:8)=IRy_med(:,2,1:3);

irf_lb(:,1)=IRg_lb(:,1,2)/max(IRg_med(:,1,2));
irf_lb(:,2:4)=IRy_lb(:,1,1:3);
irf_lb(:,5)=IRg_lb(:,2,2)/max(IRg_med(:,2,2));
irf_lb(:,6:8)=IRy_lb(:,2,1:3);

irf_ub(:,1)=IRg_ub(:,1,2)/max(IRg_med(:,1,2));
irf_ub(:,2:4)=IRy_ub(:,1,1:3);
irf_ub(:,5)=IRg_ub(:,2,2)/max(IRg_med(:,2,2));
irf_ub(:,6:8)=IRy_ub(:,2,1:3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure,

x=1:lags;
for index=1:8 %negative shock first
    if index>4
        cl='red';
    else
        cl='blue';
    end
    
    subplot(2,4,index),hold on,
    
    
    error1_pos = (irf_med(:,index) - irf_lb(:,index));
    error2_pos = (irf_ub(:,index)-irf_med(:,index));
    
    g2 = plot(x,-.01+0*irf_med(:,index) ,'LineWidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    H1= shadedErrorBar(x,irf_med(:,index),[error2_pos error1_pos], cl);
    plot(irf_med(:,index),'Linewidth',2,'Color',cl);
    
    xlim([1 30]);
    if index==1 |index==5
        title('G')
        ylim([-.01 1.5])
    end
    
    if index==2 |index==6
        title('Y, UR low')
        ylim([-.5 2])
    end
    
    if index==3 |index==7
        title('Y, UR average')
        ylim([-.5 2])
    end
    
    if index==4 |index==8
        title('Y, UR high')
        ylim([-.5 2])
    end
    
    
    if index==1
        ylabel('Expansionary shock','Fontsize',8)
    end
    if index==5
        ylabel('Contractionary shock','Fontsize',8)
    end
    set(gca,'FontSize',8)
    
    if kk==2
        ylim([0 1.25])
    end
    if kk==3
        ylim([-.5 1.5])
    end
    if kk==4
        ylim([-1 2.25])
    end
%     orient landscape
   
    set(gca,'FontSize',8)
    
    set(gca,'xticklabel',0:10:28)
    set(gca,'xtick',1:10:29)
    hold off
end

print('Figure7','-depsc') 
savefig('Figure7')

%moments of irfs
upper=95;
lower=5;
median_neg=prctile(IRF_neg,50,3); %first dimension is the observable, second the horizon
lower_neg=prctile(IRF_neg,lower,3);
upper_neg=prctile(IRF_neg,upper,3);

median_pos=prctile(IRF_pos,50,3);
lower_pos=prctile(IRF_pos,lower,3);
upper_pos=prctile(IRF_pos,upper,3);

mmax=squeeze(median(Mmax,1));
mmaxlow=squeeze(prctile(Mmax,lower,1));
mmaxhigh=squeeze(prctile(Mmax,upper,1));

msum=squeeze(median(Msum,1));
msumlow=squeeze(prctile(Msum,lower,1));
msumhigh=squeeze(prctile(Msum,upper,1));



