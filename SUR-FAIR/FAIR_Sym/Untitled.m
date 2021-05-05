
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Multipliers from LP regs obtained in Stata (Ramey-Zubairy code) 
% (all 90% intervals)
data_mRZ=csvread('../../../Ramey_Zubairy_replication_codes/RZ_M.csv',1,1);
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT IRFs of LP from Ramey's Stata code
%Normalized IRFs
x=(1:length(irfRZ_neg))';
% figure,
for jj=1:setup.dep_variables
    
%     %Plot IRFs to pos shocks
%     subplot(2,2,jj),
%     set(gca,'FontSize',8)
%     set(gca,'xtick',0:10:30)
%     hold on
%     
%     murm=mean(median_pos(jj,1:20)./max(median_pos(1,:)))-mean(irfRZ_pos(1:20,jj)./max(irfRZ_pos(:,1)));
%     
%     error1_pos = (irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)) - irfRZlb_pos(:,jj)./max(irfRZ_pos(:,1)));
%     error2_pos = (irfRZub_pos(:,jj)./max(irfRZ_pos(:,1)) - irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)));
%     H= shadedErrorBar(x,irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)),[error2_pos error1_pos], 'b');
%     
%     irfFAIR=median_pos(jj,:)./max(median_pos(1,:))-(jj==1)*murm;
%     
%     g1=plot(irfFAIR,'Linewidth',1,'Color','b','LineStyle','--')
%     plot(irfRZ_pos(:,jj)./max(irfRZ_pos(:,1)),'Linewidth',2,'Color','b','LineStyle','-');
%     
%     if jj==1
%         ylabel('Expansionary shock','FontSize',10)
%         title('Government spending','FontSize',10)
%         legend([g1],{'FAIR'}), legend('boxoff')
%         ylim([-.0 1.5])
%     else
%         title('Output')
%         ylim([0 1.5])
%     end
%     xlim([0.7 hor])
    
    
    %Plot IRFs to neg shocks
    subplot(3,3,6+jj), 
    set(gca,'FontSize',8)
    set(gca,'xtick',0:10:30)
    hold on
    
    murm=mean(median_neg(jj,1:20)./max(median_neg(1,:)))-mean(irfRZ(1:20,jj)./max(irfRZ(:,1)));
    
    error1_neg = (irfRZ(:,jj)./max(irfRZ(:,1)) - irfRZlb(:,jj)./max(irfRZ(:,1)));
    error2_neg = (irfRZub(:,jj)./max(irfRZ(:,1)) - irfRZ(:,jj)./max(irfRZ(:,1)));
    H= shadedErrorBar(x,irfRZ(:,jj)./max(irfRZ(:,1)),[error2_neg error1_neg], 'k');
    
    irfFAIR=median_neg(jj,:)./max(median_neg(1,:))-(jj==1)*murm;
    
    g1=plot(irfFAIR,'Linewidth',1,'Color','k','LineStyle','--')
    plot(irfRZ(:,jj)./max(irfRZ(:,1)),'Linewidth',2,'Color','k','LineStyle','-');
    
    
    if jj==1
%         ylabel('Contractionary shock','FontSize',10)
        title('Government spending','FontSize',10)
        legend([g1],{'FAIR'}), legend('boxoff')
        ylim([-.0 1.5])
        
    else
        title('Output','FontSize',10)
        ylim([-.0 1.5])
    end
    xlim([0.7 hor])
    set(gca,'xtick',0:10:30)
end

 

hh=(1:length(mRZ))';

subplot(3,3,9),
error1_neg = (mRZ - mRZlb);
error2_neg = -(mRZ - mRZub);
H= shadedErrorBar(hh,mRZ,[error2_neg  error1_neg], 'k');
hold on, plot(1+0*mRZ,'Linewidth',1,'Color',[.7 .7 .7],'LineStyle','--');
plot(MsumHor(:,1),'Linewidth',1,'Color','k','LineStyle','--');
plot(mRZ,'Linewidth',2,'Color','k');
set(gca,'FontSize',8)
ylim([0 2]), xlim([1 30]), title('$\mathcal M$'), box off
set(gca,'xtick',0:10:30)



% print('../Results/IRF_asym_Ramey_LP_1890','-depsc','-painters');
% print ('..\..\..\..\Papers\FP\Graphs\IRF_asym_Ramey_LP_1890','-depsc','-painters');
