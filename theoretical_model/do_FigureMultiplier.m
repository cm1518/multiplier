clear all; close all

lower=16;
upper=84;

load ('sumstats_Ramey2.mat','Msum','Mmax','INDICATOR')
msum_Ramey=squeeze(median(Msum,1));
INDICATOR_Ramey = INDICATOR;
msumlow_Ramey=squeeze(prctile(Msum,lower,1));
msumhigh_Ramey=squeeze(prctile(Msum,upper,1));
msumlow2_Ramey=squeeze(prctile(Msum,lower,1));
msumhigh2_Ramey=squeeze(prctile(Msum,upper,1));
error1_pos_Ramey = msum_Ramey - msumlow_Ramey;
error2_pos_Ramey = msumhigh_Ramey - msum_Ramey;


load ('sumstats_BPAG2.mat','Msum','Mmax','INDICATOR')
msum_AG=squeeze(median(Msum,3));
INDICATOR_AG = INDICATOR;
msumlow_AG=squeeze(prctile(Msum,lower,3));
msumhigh_AG=squeeze(prctile(Msum,upper,3));
msumlow2_AG=squeeze(prctile(Msum,lower,3));
msumhigh2_AG=squeeze(prctile(Msum,upper,3));
error1_pos_AG = msum_AG - msumlow_AG;
error2_pos_AG = msumhigh_AG - msum_AG;


load('sumstats_Model.mat','INDICATOR_Model','mult_exp','mult_contr')
mult_model(1,:) = mult_contr;
mult_model(2,:) = mult_exp;

%Figure for the multiplier:
figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(0,'DefaultAxesFontSize',12);
set(gcf,'units','points','position',[100,100,800,400]);

ColorOrder ={'b','r'};
pos(1,:) = [0.08 0.165 0.38 0.76];
pos(2,:) = [0.58 0.165 0.38 0.76];


pos_line = [2,1];
for j=1:2
    ha=subplot(1,2,j);
    set(ha,'position',pos(j,:));
    hold on
    
    H1 = shadedErrorBar(INDICATOR_Ramey(:),msum_Ramey(pos_line(j),:),[error2_pos_Ramey(pos_line(j),:); error1_pos_Ramey(pos_line(j),:)],[],0,0.25);    
    H2 = shadedErrorBar(INDICATOR_AG(:),msum_AG(pos_line(j),:),[error2_pos_AG(pos_line(j),:); error1_pos_AG(pos_line(j),:)],[],1,0.75);    
    
    g(1) = plot(100*INDICATOR_Model(:),mult_model(pos_line(j),:),'-o','Color','k','LineWidth',3,'MarkerSize',8);    
    g(2) = plot(INDICATOR_AG(:),msum_AG(pos_line(j),:), '-','Color',ColorOrder{j}, 'LineWidth', 2);    
    g(3) = plot(INDICATOR_Ramey(:),msum_Ramey(pos_line(j),:), '--','Color',ColorOrder{j},'LineWidth', 2);
    g(4) = plot(INDICATOR_Ramey(:),ones(length(INDICATOR_Ramey),1),'-k','LineWidth',1);
    xlim([-1.01 2]);
    xlabel('Unemployment Rate (dev. from trend)', 'Interpreter','latex','FontWeight','normal','FontSize',16)
    
    ylim([-0.01 2])
    set(gca,'xtick',-3:1:2)
    set(gca,'ytick',-1:0.25:2)
    
    if j==1
        ylabel('$\mathcal M^+$','Interpreter','latex','FontWeight','normal','FontSize',20);
        title('Expansionary shock','Interpreter','latex','FontSize',20)
    end
    
    if j==2
        ylabel('$\mathcal M^-$','Interpreter','latex','FontWeight','normal','FontSize',20);
        title('Contractionary shock','Interpreter','latex','FontSize',20)
        
        h(1) = plot(NaN,NaN,'-o','Color','k','LineWidth',3,'MarkerSize',8);
        h(2) = plot(NaN,NaN, '-','Color','k', 'LineWidth', 2);
        h(3) = plot(NaN,NaN, '--','Color','k', 'LineWidth', 2);
        
        legend1=legend(h(1:3),'Model','Data, Recursive id.','Data, Narrative id.');
        set(legend1,'Orientation','horizontal','Position',[0.18 0.008 0.7 0.06],'FontSize',14,'box','off');
    end
    hold off
end

print('Figure_10','-depsc','-painters')

savefig('Figure_10')
