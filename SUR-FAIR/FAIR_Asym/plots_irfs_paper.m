close all

setup.number_of_draws=setup.keep_draw*size(draws,2);
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,5000,1);

YGr=1; 
hor=30;

indicator=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load IRFs from Ramey-Zubairy LP regs obtained in Stata (I use their code!)
data_irfsRZ=csvread('../../Ramey_Zubairy_replication_codes/RZ_irfs.csv',1,1);
irfRZg=data_irfsRZ(:,1:3);
irfRZy=data_irfsRZ(:,4:6);
irfRZg_pos=data_irfsRZ(:,7:9);
irfRZy_pos=data_irfsRZ(:,10:12);
irfRZg_neg=data_irfsRZ(:,13:15);
irfRZy_neg=data_irfsRZ(:,16:18);

irfRZ=[irfRZg(:,1) irfRZy(:,1)];

irfRZ_pos=[irfRZg_pos(:,1) irfRZy_pos(:,1)];
irfRZub_pos=[irfRZg_pos(:,2) irfRZy_pos(:,2)];
irfRZlb_pos=[irfRZg_pos(:,3) irfRZy_pos(:,3)];

irfRZ_neg=[irfRZg_neg(:,1) irfRZy_neg(:,1)];
irfRZub_neg=[irfRZg_neg(:,2) irfRZy_neg(:,2)];
irfRZlb_neg=[irfRZg_neg(:,3) irfRZy_neg(:,3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Multipliers from LP regs obtained in Stata (Ramey-Zubairy code)
% (all 90% intervals)
data_mRZ=csvread('../../Ramey_Zubairy_replication_codes/RZ_M.csv',1,1);
%Attention, unliek with IRFs, it's ordered: pt estimate, std error
mRZ=data_mRZ(:,1);
mRZub=data_mRZ(:,1)+1.654*data_mRZ(:,2);
mRZlb=data_mRZ(:,1)-1.654*data_mRZ(:,2);

mRZ_pos=data_mRZ(:,3);
mRZ_posub=data_mRZ(:,3)+1.654*data_mRZ(:,4);
mRZ_poslb=data_mRZ(:,3)-1.654*data_mRZ(:,4);

mRZ_neg=data_mRZ(:,5);
mRZ_negub=data_mRZ(:,5)+1.654*data_mRZ(:,6);
mRZ_neglb=data_mRZ(:,5)-1.654*data_mRZ(:,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%



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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT Normalized IRFs along with cumulative multiplier:

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

x=1:setup.lags_shocks;
figure,
for jj=1:setup.dep_variables
    
    subplot(2,3,jj),
    set(gca,'FontSize',8)
    hold on
    %     g2 = plot(x,0*median_pos(2,:) ,'LineWidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    error1_pos = (median_pos(jj,:)./max(median_pos(1,:)) - lower_pos(jj,:)./max(median_pos(1,:)))';
    error2_pos = (upper_pos(jj,:)./max(median_pos(1,:)) - median_pos(jj,:)./max(median_pos(1,:)))';
    H= shadedErrorBar(x,median_pos(jj,:)./max(median_pos(1,:)),[error2_pos'; error1_pos'], 'b');
    g1=plot(irfsym(jj,:)./max(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--');
    %     murm=mean(median_pos(jj,1:hor)./max(median_pos(1,:)))-mean(irfRZ_pos(1:hor,jj)./max(irfRZ_pos(:,1)));
    %     g1=plot(irfRZ_pos(:,jj)./max(irfRZ_pos(:,1))+(jj==1)*murm,'Linewidth',1,'Color','b','LineStyle','--');
    plot(median_pos(jj,:)./max(median_pos(1,:)),'Linewidth',2,'Color','b')
    
    if jj==1
        ylabel('Expansionary shock','Fontsize',10)
        title('Government spending','Fontsize',10)
        legend([g1],{'Linear'},'FontSize',8), legend('boxoff')
        ylim([-.0 1.5])
    else
        title('Output','FontSize',10)
        ylim([0 1.5])
    end
    xlim([1 30]), set(gca,'xtick',0:10:30)
    
    subplot(2,3,3+jj),
    set(gca,'FontSize',8)
    hold on
    %     g2 = plot(x,0*median_neg(2,:) ,'LineWidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    error1_neg = (median_neg(jj,:)./max(median_neg(1,:)) - lower_neg(jj,:)./max(median_neg(1,:)))';
    error2_neg = (upper_neg(jj,:)./max(median_neg(1,:)) - median_neg(jj,:)./max(median_neg(1,:)))';
    H= shadedErrorBar(x,median_neg(jj,:)./max(median_neg(1,:)),[error2_neg'; error1_neg'], 'r');
    g1=plot(irfsym(jj,:)./max(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--');
    %     murm=mean(median_neg(jj,1:hor)./max(median_neg(1,:)))-mean(irfRZ_neg(1:hor,jj)./max(irfRZ_neg(:,1)));
    %     g1=plot(irfRZ_neg(:,jj)./max(irfRZ_neg(:,1))-(jj==1)*murm,'Linewidth',1,'Color','r','LineStyle','--');
    plot(median_neg(jj,:)./max(median_neg(1,:)),'Linewidth',2,'Color','r')
    
    %     plot(1:setup.lags_shocks,median_neg(jj,:)./max(median_neg(1,:)),'Linewidth',2,'Color','r')
    %     hold on, plot(1:setup.lags_shocks,lower_neg(jj,:)./max(median_neg(1,:)),1:setup.lags_shocks,upper_neg(jj,:)./max(median_neg(1,:)),'Color','r')
    %     plot(irfsym(jj,:)./max(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--')
    
    if jj==1
        ylabel('Contractionary shock','FontSize',10)
        title('(-) Government spending','FontSize',10)
        legend([g1],{'Linear'},'FontSize',8), legend('boxoff')
        ylim([-.0 1.5])
        
             
    else
        title('(-) Output','FontSize',10)
        ylim([-.0 1.5])
    end
    xlim([1 30]), set(gca,'xtick',0:10:30)
    
end

%plot cumulative multiplier by horizon
hh=(1:hor)';

subplot(2,3,3),
error1_pos = (MsumHor(:,2) - MsumHor_lb(:,2));
error2_pos = -(MsumHor(:,2) - MsumHor_ub(:,2));
H= shadedErrorBar(hh,MsumHor(:,2),[error2_pos  error1_pos], 'b');
hold on, plot(1+0*MsumHor(:,1),'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
plot(cumsum(irfsym(2,:))./cumsum(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--')
plot(MsumHor(:,2),'Linewidth',2,'Color','b');
set(gca,'FontSize',8)
ylim([0 2]), xlim([1 30]), title('$\mathcal M^+$'), box off
set(gca,'xtick',0:10:30)

subplot(2,3,6),
error1_neg = (MsumHor(:,1) - MsumHor_lb(:,1));
error2_neg = -(MsumHor(:,1) - MsumHor_ub(:,1));
H= shadedErrorBar(hh,MsumHor(:,1),[error2_neg  error1_neg], 'r');
hold on, plot(1+0*MsumHor(:,1),'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
plot(cumsum(irfsym(2,:))./cumsum(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--')
plot(MsumHor(:,1),'Linewidth',2,'Color','r');
set(gca,'FontSize',8)
ylim([0 2]), xlim([1 30]), title('$\mathcal M^-$'), box off
set(gca,'xtick',0:10:30)

print('Figure3','-depsc') 
savefig('Figure3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT IRFs of LP from Ramey's Stata code:3*3 figure with symmetric IRFs
%Normalized IRFs
x=(1:length(irfRZ_neg))';
fig1=openfig('IRFs_LP_FAIR.fig');

for jj=1:setup.dep_variables
    
    %Plot IRFs to pos shocks
    subplot(3,3,3+jj),
    set(gca,'FontSize',8)
    set(gca,'xtick',0:10:30)
    hold on, plot(0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    
    
    murm=-.2+mean(median_pos(jj,1:hor)./max(median_pos(1,:)))-mean(irfRZ_pos(1:hor,jj)./max(irfRZ_pos(:,1)));
    
    error1_pos = (irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)) - irfRZlb_pos(:,jj)./max(irfRZ_pos(:,1)));
    error2_pos = (irfRZub_pos(:,jj)./max(irfRZ_pos(:,1)) - irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)));
    H= shadedErrorBar(x,irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)),[error2_pos error1_pos], 'b');
    
    irfFAIR=median_pos(jj,:)./max(median_pos(1,:))-(jj==1)*murm;
    
    g1=plot(irfFAIR,'Linewidth',1,'Color','b','LineStyle','--');
    plot(irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)),'Linewidth',2,'Color','b','LineStyle','-');
    
    if jj==1
        ylabel('Exp. shock','FontSize',10)
        title('G','FontSize',10)
        legend([g1],{'FAIR'}), legend('boxoff')
        ylim([-.5 1.5])
    else
        title('Y','FontSize',10)
        ylim([-.5 1.5])
    end
    xlim([0.7 hor])
        set(gca,'xtick',0:10:30)
    
    
    %Plot IRFs to neg shocks
    subplot(3,3,6+jj),
    set(gca,'FontSize',8)
    set(gca,'xtick',0:10:30)
    hold on, plot(0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    
    murm=mean(median_neg(jj,1:hor)./max(median_neg(1,:)))-mean(irfRZ_neg(1:hor,jj)./max(irfRZ_neg(:,1)));
    
    error1_neg = (irfRZ_neg(:,jj)./max(irfRZ_neg(:,1)) - irfRZlb_neg(:,jj)./max(irfRZ_neg(:,1)));
    error2_neg = (irfRZub_neg(:,jj)./max(irfRZ_neg(:,1)) - irfRZ_neg(:,jj)./max(irfRZ_neg(:,1)));
    H= shadedErrorBar(x,irfRZ_neg(:,jj)./max(irfRZ_neg(:,1)),[error2_neg error1_neg], 'r');
    
    irfFAIR=median_neg(jj,:)./max(median_neg(1,:))-(jj==1)*murm;
    
    g1=plot(irfFAIR,'Linewidth',1,'Color','r','LineStyle','--');
    plot(irfRZ_neg(:,jj)./max(irfRZ_neg(:,1)),'Linewidth',2,'Color','r','LineStyle','-');
    
    
    if jj==1
        ylabel('Contr. shock','FontSize',10)
        title('(-) G','FontSize',10)
        legend([g1],{'FAIR'}), legend('boxoff')
        ylim([-.5 1.5])
        
    else
        title('(-) Y','FontSize',10)
        ylim([-.5 1.5])
    end
    xlim([0.7 hor])
    set(gca,'xtick',0:10:30)
end


hh=(1:length(mRZ))';

subplot(3,3,6),
error1_pos = (mRZ_pos - mRZ_poslb);
error2_pos = -(mRZ_poslb - mRZ_pos);
H= shadedErrorBar(hh,mRZ_pos,[error2_pos  error1_pos], 'b');
hold on, plot(1+0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
plot(MsumHor(:,2),'Linewidth',1,'Color','b','LineStyle','--');
plot(mRZ_pos,'Linewidth',2,'Color','b');
set(gca,'FontSize',8)
ylim([0 2]), xlim([1 24]), title('$\mathcal M^+$'), box off
set(gca,'xticklabel',0:10:30,'xtick',1:11:31)



subplot(3,3,9),
error1_neg = (mRZ_neg - mRZ_neglb);
error2_neg = -(mRZ_neglb - mRZ_neg);
H= shadedErrorBar(hh,mRZ_neg,[error2_neg  error1_neg], 'r');
hold on, plot(1+0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
plot(MsumHor(:,1),'Linewidth',1,'Color','r','LineStyle','--');
plot(mRZ_neg,'Linewidth',2,'Color','r');
set(gca,'FontSize',8)
ylim([0 2]), xlim([1 24]), title('$\mathcal M^-$'), box off
set(gca,'xticklabel',0:10:30,'xtick',1:11:31)


print('Figure4','-depsc') 
savefig('Figure4')

 



