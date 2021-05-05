

setup.number_of_draws=setup.keep_draw*length(draws);

inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,1000,1);

variablename={'Government spending','Tax','Output'};
% Horizon of IRF
setup.irfhor=30;

YGr=1;

%Error bands confidence:
upper=95;
lower=5;

%Horizon for msum
H=setup.irfhor;
H=20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
    index_ind=0;
    for indicator=[-1 0 2];
        
        index_ind=index_ind+1;
        
        % if exist('IRFdist.mat')~=2
        
        MsumS=zeros(2,length(inds));MmaxS=zeros(2,length(inds));
        
        for kk=1:length(indices_for_draws)
            
            index=0;
            for size_of_shock=1:-2:-1;
                index=index+1;
                shock_contemp=zeros(setup.size_obs,1);
                shock_contemp(setup.index_unrestricted)=size_of_shock;
                ind_vec=indicator*[ones(1,setup.lags+1)];
                for jj=1:setup.lags+1
                    epsilon=[zeros(setup.size_obs,jj-1) shock_contemp zeros(setup.size_obs,setup.lags+1-jj)];
                    [ Sigma, intercept] = unwrap_NL_IRF( draws(:,indices_for_draws(kk)),epsilon,setup ,ind_vec,0);
                    if jj==1
                        ind_vec=[ind_vec(1,1:end-1) indicator];
                    else
                        ind_vec=[ind_vec(1,2:end) 0];
                    end
                    
                    IRFs(:,jj,kk,index)=Sigma(:,setup.index_unrestricted,jj);
                    
                end
                
                %Re-scale by size of G shock
                IRFs([1 3 4],:,kk,index)=IRFs([1 3 4],:,kk,index)/max(abs(IRFs(2,:,kk,index)));
                %IRFs(:,:,kk,index)=IRFs(:,:,kk,index)/max(abs(IRFs(2,:,kk,index)));
                
                clear ind_vec
            end
            
            
        end
        
        IRFsmed=squeeze(prctile(IRFs,50,3));
        IRFslow=squeeze(prctile(IRFs,lower,3));
        IRFshigh=squeeze(prctile(IRFs,upper,3));
        
        
        IRy_med(:,:,index_ind)=squeeze(IRFsmed(4,:,:));
        IRy_lb(:,:,index_ind)=squeeze(IRFslow(4,:,:));
        IRy_ub(:,:,index_ind)=squeeze(IRFshigh(4,:,:));
        
        
        IRg_med(:,:,index_ind)=squeeze(IRFsmed(2,:,:));
        IRg_lb(:,:,index_ind)=squeeze(IRFslow(2,:,:));
        IRg_ub(:,:,index_ind)=squeeze(IRFshigh(2,:,:));
        
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

    save irfs_state.mat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure,

x=1:setup.lags+1;
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
        ylim([-.01 1.2])
    end
    
    if index==2 |index==6
        title('Y, UR low')
        ylim([-.5 3])
    end
    
    if index==3 |index==7
        title('Y, UR average')
        ylim([-.5 3])
    end
    
    if index==4 |index==8
        title('Y, UR high')
        ylim([-.5 3])
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
   
    set(gca,'xticklabel',0:10:28)
    set(gca,'xtick',1:10:29)
    hold off
end



print('Figure6','-depsc') 
savefig('Figure6')


