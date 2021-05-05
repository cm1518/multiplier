% close all
% clear all
% load ('sumstats.mat','Msum','Mmax','INDICATOR')

% std_ind=std(setup.indicator);
% mu_ind=mean(setup.indicator);
%
% % INDICATOR=mu_ind+(-1.5:.5:2)*std_ind;
% INDICATOR=-1:.5:2;

msum=squeeze(median(Msum,1))';
msumlow=squeeze(prctile(Msum,lower,1))';
msumhigh=squeeze(prctile(Msum,upper,1))';
msumlow2=squeeze(prctile(Msum,lower2,1))';
msumhigh2=squeeze(prctile(Msum,upper2,1))';

mmax=squeeze(median(Mmax,1))';
mmaxlow=squeeze(prctile(Mmax,lower,1))';
mmaxhigh=squeeze(prctile(Mmax,upper,1))';
mmaxlow2=squeeze(prctile(Mmax,lower2,1))';
mmaxhigh2=squeeze(prctile(Mmax,upper2,1))';

msum

mmax

% Linear multiplier:
% TBW

% 
% %% FIGURE with shaded error-bands
% 
% % Mmax
% hFig = figure(1);
% % set(hFig, 'Position', [300 200 1000 500]) %x,y, with, height
% %Figure for error distribution
% for j=1:2
%     if j==1,
%         h3=subplot(223); hold on;
%         %         plot(xout,[.5*Ap*exp(-(xout-mup).^2/(2*sigmap^2)); ],'LineWidth',1.5)
%         bar(xout,.5*r_BMpos_binv./sum(r_BMpos_binv),1) %normailzing the number of shocks to 1
% %              bar(xout,.5*r_BMpos_binv,1) %NOT normailzing the number of shocks to 1
%         xlim([-1.5 2.5]);
%         ylim([0 .1])
% %         ylim([0 50])
%         %   xlabel('Unemployment rate')
%         %         ylim([-.1 .1])
%         xlabel('(detrended) UR')
%         ax=get(h3,'Position');
%         ax(2)=ax(2)-.0;
%         ax(4)=ax(4)-0.20; %or wathever
%         set(h3,'Position',ax);
%         set(gca,'ytick',[]), set(gca,'Color','none');
%         %         set(gca,'ycolor',get(gcf,'color'));
%         set(gca,'xtick',-3:1:3)
%         ylabel('  Shocks\newlineDistribution','Fontsize',8)
%     end
%     if j==2
%         h=subplot(224); hold on;
%         %         plot(xout,[.5*Ap*exp(-(xout-mup).^2/(2*sigmap^2)); ],'LineWidth',1.5)
%         bar(xout,.5*r_BMneg_binv./sum(r_BMneg_binv),1) %normailzing the number of shocks to 1
% %               bar(xout,.5*r_BMneg_binv,1) %NOT normailzing the number of shocks to 1
%         xlim([-1.5 2.5]);
%         ylim([0 .1])
% %         ylim([0 50])
%         xlabel('(detrended) UR')
%         ax=get(h,'Position');
%         ax(2)=ax(2)-.0;
%         ax(4)=ax(4)-0.20; %or wathever
%         set(h,'Position',ax);
%         set(gca,'ytick',[]), set(gca,'Color','none');
%         %         set(gca','Xticklabel',-1.5:1:3.5)
%         set(gca,'xtick',-1:1:3)
%         %         set(gca,'ycolor',get(gcf,'color'));
%         %         ylabel('  Shocks\newlineDistribution','Fontsize',8)
%     end
% end
% 
% 
% 
% %Plot mmax
% mu_ur = 0;
% x = mu_ur+INDICATOR;
% for j=1:2
%     error1_pos = (mmax(:,j) - mmaxlow(:,j))';
%     error2_pos = (mmaxhigh(:,j) - mmax(:,j))';
%     
%     error1_pos2 = (mmax(:,j) - mmaxlow2(:,j))';
%     error2_pos2 = (mmaxhigh2(:,j) - mmax(:,j))';
% 
%     if j==1,
%         h1=subplot(222); hold on
%         H2= shadedErrorBar_RB(x,mmax(:,j),[error2_pos2; error1_pos2], 'r',0,.5);
%         H= shadedErrorBar_RB(x,mmax(:,j),[error2_pos; error1_pos],'r',0,1);
%         g1 = plot(x,mmax(:,j), 'r', 'LineWidth', 2);
%         g2 = plot(x,1+0*mmax(:,1),'--k','LineWidth',2);
%         %         xlim([min(x)*.995 max(x)*1.0]);
%         xlim([-1.5 2.5]);
% %         legend([g1],{'FP'},'Location','NorthWest'), legend('boxoff')
%         title('Contractionary shock')
%         ylim([0 3])
%         ax=get(h1,'Position');
%         ax(2)=ax(2)-.17;
%         ax(4)=ax(4)+0.18; %or wathever
%         set(h1,'Position',ax);
%         set(gca,'xtick',-3:1:3)
%         set(gca,'ytick',-1:1:3)
%     end
%     if j==2
%         h2=subplot(221); hold on
%         H2= shadedErrorBar_RB(x,mmax(:,j),[error2_pos2; error1_pos2], 'b',0,.5);
%         H= shadedErrorBar_RB(x,mmax(:,j),[error2_pos; error1_pos],'b',0,1);
%         g1 = plot(x,mmax(:,j), 'b', 'LineWidth', 2);
%         g2 = plot(x,1+0*mmax(:,1),'--k','LineWidth',2);
%         %         xlim([min(x)*.995 max(x)*1.0]);
%         xlim([-1.5 2.5]);
%         ylim([0 3])
%         ylabel('Multiplier (max)')
%         legend([g1],{'FP'},'Location','NorthWest'), legend('boxoff')
%         title('Expansionary shock')
%         ax=get(h2,'Position');
%         ax(2)=ax(2)-.17;
%         ax(4)=ax(4)+0.18; %or wathever
%         set(h2,'Position',ax);
%         set(gca,'xtick',-1:1:3)
%         set(gca,'ytick',-1:1:3)
%     end
% end
% 
% 
% 
% 
% 
% %% Msum
% hFig = figure(2);
% % set(hFig, 'Position', [300 200 1000 500]) %x,y, with, height
% %Figure for error distribution
% for j=1:2
%     if j==1,
%         h3=subplot(223); hold on;
%         %         plot(xout,[.5*Ap*exp(-(xout-mup).^2/(2*sigmap^2)); ],'LineWidth',1.5)
%         bar(xout,.5*r_BMpos_binv./sum(r_BMpos_binv),1)
%         xlim([-1.5 2.5]);
%         ylim([0 .1])
%         %   xlabel('Unemployment rate')
%         %         ylim([-.1 .1])
%         xlabel('(detrended) UR')
%         ax=get(h3,'Position');
%         ax(2)=ax(2)-.0;
%         ax(4)=ax(4)-0.20; %or wathever
%         set(h3,'Position',ax);
%         set(gca,'ytick',[]), set(gca,'Color','none');
%         %         set(gca,'ycolor',get(gcf,'color'));
%         set(gca,'xtick',-3:1:3)
%         ylabel('  Shocks\newlineDistribution','Fontsize',8)
%     end
%     if j==2
%         h=subplot(224); hold on;
%         %         plot(xout,[.5*Ap*exp(-(xout-mup).^2/(2*sigmap^2)); ],'LineWidth',1.5)
%         bar(xout,.5*r_BMneg_binv./sum(r_BMneg_binv),1)
%         xlim([-1.5 2.5]);
%         ylim([0 .1])
%         xlabel('(detrended) UR')
%         ax=get(h,'Position');
%         ax(2)=ax(2)-.0;
%         ax(4)=ax(4)-0.20; %or wathever
%         set(h,'Position',ax);
%         set(gca,'ytick',[]), set(gca,'Color','none');
%         %         set(gca','Xticklabel',-1.5:1:3.5)
%         set(gca,'xtick',-1:1:3)
%         %         set(gca,'ycolor',get(gcf,'color'));
%         %         ylabel('  Shocks\newlineDistribution','Fontsize',8)
%     end
% end
% 
% %Plot msum
% mu_ur = 0;
% x = mu_ur+INDICATOR;
% for j=1:2
%     error1_pos = (msum(:,j) - msumlow(:,j))';
%     error2_pos = (msumhigh(:,j) - msum(:,j))';
%     
%     error1_pos2 = (msum(:,j) - msumlow2(:,j))';
%     error2_pos2 = (msumhigh2(:,j) - msum(:,j))';
%      
%     if j==1,
%         h1=subplot(222); hold on;
%         H2= shadedErrorBar_RB(x,msum(:,j),[error2_pos2; error1_pos2], 'r',0,.5);
%         H= shadedErrorBar_RB(x,msum(:,j),[error2_pos; error1_pos],'r',0,1);
%         g1 = plot(x,msum(:,j), 'r', 'LineWidth', 2);
%         g2 = plot(x,1+0*msum(:,1),'--k','LineWidth',2);
%         %         xlim([min(x)*.995 sum(x)*1.0]);
%         xlim([-1.5 2.5]);
% %         legend([g1],{'FP'},'Location','NorthWest'), legend('boxoff')
%         title('Contractionary shock')
%         ylim([0 3])
%         ax=get(h1,'Position');
%         ax(2)=ax(2)-.17;
%         ax(4)=ax(4)+0.18; %or wathever
%         set(h1,'Position',ax);
%         set(gca,'xtick',-3:1:3)
%         set(gca,'ytick',-1:1:3)
%     end
%     if j==2
%         h2=subplot(221); hold on;
%          H2= shadedErrorBar_RB(x,msum(:,j),[error2_pos2; error1_pos2], 'b',0,.5);
%         H= shadedErrorBar_RB(x,msum(:,j),[error2_pos; error1_pos],'b',0,1);
%         g1 = plot(x,msum(:,j), 'b', 'LineWidth', 2);
%         g2 = plot(x,1+0*msum(:,1),'--k','LineWidth',2);
%         %         xlim([min(x)*.995 sum(x)*1.0]);
%         xlim([-1.5 2.5]);
%         ylim([0 3]);
%         ylabel('Multiplier (sum)')
% %         legend([g1],{'FP'},'Location','NorthWest'), legend('boxoff')
%         title('Expansionary shock')
%         ax=get(h2,'Position');
%         ax(2)=ax(2)-.17;
%         ax(4)=ax(4)+0.18; %or wathever
%         set(h2,'Position',ax);
%         set(gca,'xtick',-1:1:3)
%         set(gca,'ytick',-1:1:3)
%     end
% end
% 
% 
% 
