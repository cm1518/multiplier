


load results_asym


load seed_for_fig2 
rng(s);



setup.number_of_draws=setup.keep_draw*length(draws);

inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,1000,1);

set(0,'defaulttextInterpreter','latex')

setup.irfhor=30;

YGr=1;

%Error bands confidence:
upper=95;
lower=5;

%Horizon for msum
H=setup.irfhor;
H=20;

%Load results from symmetric (linear) FAIR:
load IRsym.mat Msum_fairsym IR_fairsym

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indicator=0;

    
    MsumS=zeros(2,length(inds));MmaxS=zeros(2,length(inds));
    
    for kk=1:length(indices_for_draws)
        
        index=0;
        for size_of_shock=1:-2:-1;
            index=index+1;
            shock_contemp=zeros(setup.size_obs,1);
            shock_contemp(setup.index_unrestricted)=size_of_shock;
            ind_vec=[zeros(1,setup.lags)];
            for jj=1:setup.lags+1
                epsilon=[zeros(setup.size_obs,jj-1) shock_contemp zeros(setup.size_obs,setup.lags+1-jj)];
                [ Sigma, intercept] = unwrap_NL_IRF( draws(:,indices_for_draws(kk)),epsilon,setup ,ind_vec);
                if jj==1
                    ind_vec=[ind_vec(1,1:end-1) indicator];
                else
                    ind_vec=[ind_vec(1,2:end) 0];
                end
                
                IRFs(:,jj,kk,index)=Sigma(:,setup.index_unrestricted,jj);
                
            end
            for h=1:H
                MsumS(index,kk,h)=sum(IRFs(4,1:h+1,kk,index),2)./sum(IRFs(2,1:h+1,kk,index),2);
            end
            MmaxS(index,kk)=max(IRFs(4,:,kk,index))./max(IRFs(2,:,kk,index));
            clear ind_vec
        end
        
        
    end
    
    IRFsmed=squeeze(median(IRFs,3));
%     IRFsmed=squeeze(mean(IRFs,3));
    IRFslow=squeeze(prctile(IRFs,lower,3));
    IRFshigh=squeeze(prctile(IRFs,upper,3));
    
    save ('IRFdist.mat', 'IRFs', 'IRFsmed',  'IRFslow', 'IRFshigh','MsumS','MmaxS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure for paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


counter=1;


for kk=1:setup.horizon+1
    setup.store_responses(:,:,kk)=VAR_MA_implied_higher_order( setup,kk-1);
end
store_responses=setup.store_responses;


resp_1_1_VAR=store_responses(1,1,:);
resp_1_2_VAR=store_responses(1,2,:);
resp_1_3_VAR=store_responses(1,3,:);
resp_1_4_VAR=store_responses(1,4,:);
resp_2_1_VAR=store_responses(2,1,:);
resp_2_2_VAR=store_responses(2,2,:);
resp_2_3_VAR=store_responses(2,3,:);
resp_2_4_VAR=store_responses(2,4,:);
resp_3_1_VAR=store_responses(3,1,:);
resp_3_2_VAR=store_responses(3,2,:);
resp_3_3_VAR=store_responses(3,3,:);
resp_3_4_VAR=store_responses(3,4,:);
resp_4_1_VAR=store_responses(4,1,:);
resp_4_2_VAR=store_responses(4,2,:);
resp_4_3_VAR=store_responses(4,3,:);
resp_4_4_VAR=store_responses(4,4,:);

irfG0=max(resp_2_2_VAR);





variablename={'Government spending','Tax','Output'};

%moments of irfs
upper=95;
lower=5;

 
MsumHor=squeeze(prctile(MsumS,50,2))';
MsumHor_ub=squeeze(prctile(MsumS,upper,2))';
MsumHor_lb=squeeze(prctile(MsumS,lower,2))';


MsumHor_VAR=cumsum(squeeze(resp_4_2_VAR(:,:,:)))./cumsum(squeeze(resp_2_2_VAR(:,:,:)));



x=1:setup.lags;

% Plot 2*2 figure (G and Y only, like Ramey) with Cumulative Multiplier
figure,
counter=1;
x=1:setup.lags+1;
for index=1:2 %negative shock first
    if index==2
        cl='red';
    else
        cl='blue';
    end
    for kk=[2 4] %I do not plot the IRF of Tax or GE (SPF growth forecast of G)
    if counter==3
        counter=counter+1;
    end
    
        subplot(2,3,counter), hold on,
       
        if kk==2
            scl_ub=max(squeeze(IRFsmed(2,:,index)));
            scl_lb=max(squeeze(IRFsmed(2,:,index)));
            scl_m=max(squeeze(IRFsmed(2,:,index)));
        end
       
        if kk==4
            scl_ub=scl_ub/YGr;
            scl_lb=scl_lb/YGr;
            scl_m=scl_m/YGr;
        end
        
        error1_pos = (squeeze(IRFsmed(kk,:,index))/scl_m - squeeze(IRFslow(kk,:,index))/scl_lb)';
        error2_pos = (squeeze(IRFshigh(kk,:,index))/scl_ub - squeeze(IRFsmed(kk,:,index))/scl_m)';
        
        g2 = plot(x,-.01+0*squeeze(IRFsmed(kk,:,index))/scl_m ,'LineWidth',1,'Color',[.7 .7 .7],'LineStyle','--');
        H1= shadedErrorBar(x,squeeze(IRFsmed(kk,:,index))/scl_m,[error2_pos'; error1_pos'], cl);
        plot(squeeze(IRFsmed(kk,:,index))/scl_m,'Linewidth',2,'Color',cl)

        g1=plot(IR_fairsym(:,kk),'Linewidth',2,'Color','k','LineStyle','--');
        
        xlim([1 20]);
        if index==2
            str=sprintf(['(-) ',char(variablename(kk-1))]);
        else
            str=sprintf(char(variablename(kk-1)));
            
        end
        if counter==1
            ylabel('Expansionary shock','Fontsize',10)
            legend([g1],{'Linear'},'FontSize',8,'Location','SouthWest')
            legend('boxoff')
        end
        
        if counter==4
            legend([g1],{'Linear'},'FontSize',8,'Location','SouthWest')
            legend('boxoff')
            ylabel('Contractionary shock','Fontsize',10)
        end
        set(gca,'FontSize',8)
        
        if kk==2
            ylim([0 1.25])
        end
        if kk==3
            ylim([-.5 1.5])
        end
        if kk==4
            ylim([-.5 2.25])
        end
        
        title(str,'Fontsize',10)
        
        hold off
        counter=counter+1;
    end
    
end

%plot cumulative multiplier by horizon
hh=(1:H)';

subplot(2,3,3),
error1_pos = (MsumHor(:,1) - MsumHor_lb(:,1));
error2_pos = -(MsumHor(:,1) - MsumHor_ub(:,1));
shadedErrorBar(hh,MsumHor(:,1),[error2_pos  error1_pos], 'b');
hold on, plot(1+0*MsumHor(:,1),'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');

plot(Msum_fairsym(:,1),'Linewidth',2,'Color','k','LineStyle','--')
plot(MsumHor(:,1),'Linewidth',2,'Color','b');
set(gca,'FontSize',8)
ylim([0 3]), xlim([1 20]), title('$\mathcal M^+$'), box off


subplot(2,3,6),
error1_neg = (MsumHor(:,2) - MsumHor_lb(:,2));
error2_neg = -(MsumHor(:,2) - MsumHor_ub(:,2));
shadedErrorBar(hh,MsumHor(:,2),[error2_neg  error1_neg], 'r');
hold on, plot(1+0*MsumHor(:,1),'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');

plot(Msum_fairsym(:,1),'Linewidth',2,'Color','k','LineStyle','--')
plot(MsumHor(:,2),'Linewidth',2,'Color','r');
set(gca,'FontSize',8)
ylim([0 3]),  xlim([1 20]), title('$\mathcal M^-$'), box off
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display multipliers:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m_sum:
msum=prctile(MsumS,50,2);
msum_ub=(prctile(MsumS,upper,2));
msum_lb=(prctile(MsumS,lower,2));








print('Figure2','-depsc') 
savefig('Figure2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test that m_neg>m_pos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Test_statedep_fig_2

