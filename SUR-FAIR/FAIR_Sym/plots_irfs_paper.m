

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load IRFs from LP regs obtained in Stata (all 90% intervals)
data_irfsRZ=csvread('../../Ramey_Zubairy_replication_codes/RZ_irfs.csv',1,1);
irfRZg=data_irfsRZ(:,1:3);
irfRZy=data_irfsRZ(:,4:6);
irfRZg_pos=data_irfsRZ(:,7:9);
irfRZy_pos=data_irfsRZ(:,10:12);
irfRZg_neg=data_irfsRZ(:,13:15);
irfRZy_neg=data_irfsRZ(:,16:18);

irfRZ=[irfRZg(:,1) irfRZy(:,1)];
irfRZub=[irfRZg(:,2) irfRZy(:,2)];
irfRZlb=[irfRZg(:,3) irfRZy(:,3)];

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

load SUR_sym_paper;

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

MsumHor=[];
MsumHor=squeeze(msum)';
MsumHor_lb=squeeze(msum_lb)';
MsumHor_ub=squeeze(msum_ub)';

x=1:setup.lags_shocks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT IRFs of LP from Ramey's Stata code
%Normalized IRFs
x=(1:length(irfRZ_neg))';
fig1=figure,
for jj=1:setup.dep_variables
    
    %Plot IRFs to neg (symmetric) shocks
    subplot(3,3,jj), 
    set(gca,'FontSize',8)
    set(gca,'xtick',0:10:30)
    hold on, plot(0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
      
    murm=mean(median_neg(jj,1:20)./max(median_neg(1,:)))-mean(irfRZ(1:20,jj)./max(irfRZ(:,1)));
    
    error1_neg = (irfRZ(:,jj)./max(irfRZ(:,1)) - irfRZlb(:,jj)./max(irfRZ(:,1)));
    error2_neg = (irfRZub(:,jj)./max(irfRZ(:,1)) - irfRZ(:,jj)./max(irfRZ(:,1)));
    H= shadedErrorBar(x,irfRZ(:,jj)./max(irfRZ(:,1)),[error2_neg error1_neg], 'k');
    
    irfFAIR=median_neg(jj,:)./max(median_neg(1,:))-(jj==1)*murm;
    
    g1=plot(irfFAIR,'Linewidth',1,'Color','k','LineStyle','--')
    plot(irfRZ(:,jj)./max(irfRZ(:,1)),'Linewidth',2,'Color','k','LineStyle','-');
    
    
    if jj==1
%         ylabel('Contractionary shock','FontSize',10)
        title('G','FontSize',10)
        legend([g1],{'FAIR'},'FontSize',8), legend('boxoff')
        ylim([-.5 1.5])
        ylabel('Linear','FontSize',10)
        
    else
        title('Y','FontSize',10)
        ylim([-.5 1.5])
    end
    xlim([0.7 hor])
    set(gca,'xtick',0:10:30)
end

 

hh=(1:length(mRZ))';

subplot(3,3,3),
error1_neg = (mRZ - mRZlb);
error2_neg = -(mRZ - mRZub);
H= shadedErrorBar(hh,mRZ,[error2_neg  error1_neg], 'k');
hold on, plot(1+0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
plot(MsumHor(:,1),'Linewidth',1,'Color','k','LineStyle','--');
plot(mRZ,'Linewidth',2,'Color','k');
set(gca,'FontSize',8)
ylim([0 2]), xlim([1 24]), title('$\mathcal M$'), box off
set(gca,'xticklabel',0:10:30,'xtick',1:11:31)


% Save to use later in the asym case
saveas(fig1,'../FAIR_Asym/IRFs_LP_FAIR.fig');


