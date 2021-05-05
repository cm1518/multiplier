close all
clear all
load  test_NL_112615 %for detrended UR as indicator
load ('sumstats_112615.mat','Msum','Mmax','INDICATOR')
% load  test_NL_020116_Ygap %for Ygap as indicator
% load ('sumstats_Yap.mat','Msum','Mmax','INDICATOR')

% 
% std_ind=std(setup.indicator);
% mu_ind=mean(setup.indicator);
% 
% % INDICATOR=mu_ind+(-1.5:.5:2)*std_ind;
% INDICATOR=-1:.5:2;

lower=15;
upper=85;

YGr=5;

msum=squeeze(median(Msum,3));
msumlow=squeeze(prctile(Msum,lower,3));
msumhigh=squeeze(prctile(Msum,upper,3));

mmax=squeeze(median(Mmax,3));
mmaxlow=squeeze(prctile(Mmax,lower,3));
mmaxhigh=squeeze(prctile(Mmax,upper,3));

msum

mmax

%VAR-based multipliers:
for jj=1:setup.size_obs
    irf_var(:,jj)=squeeze(setup.store_responses(jj,setup.index_unrestricted,:));
end
msum_var=sum(irf_var(1:21,4))./sum(irf_var(1:21,2))*YGr;
mmax_var=max(irf_var(1:21,4))./max(irf_var(1:21,2))*YGr;



%% FIGURE with shaded error-bands

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMAX:
mu_ur = 5.6;
x = mu_ur+INDICATOR;
figure,
for j=1:2
    error1_pos = mmax(j,:) - mmaxlow(j,:);
    error2_pos = mmaxhigh(j,:) - mmax(j,:);
    if j==1,
        subplot(212), hold on
        H= shadedErrorBar(x,mmax(j,:),[error2_pos; error1_pos], 'r');
        g1 = plot(x,mmax(j,:), 'r', 'LineWidth', 2);
        g2 = plot(x,mmax_var+0*mmax(1,:),'--k','LineWidth',2);
%         xlim([min(x)*.995 max(x)*1.0]);
        xlim([5 max(x)*1.0]);
        legend([g2 g1],{'VAR','GMA'},'Location','NorthWest'), legend('boxoff')     
        ylabel('Consolidation shock')
        ylim([0 2.5])
        xlabel('Unemployment rate')
        set(gca,'ytick',0:.5:2.5)
    end
    if j==2
        subplot(211), hold on
        H= shadedErrorBar(x,mmax(j,:),[error2_pos; error1_pos], 'b')
        g1 = plot(x,mmax(j,:), 'b', 'LineWidth', 2);
        g2 = plot(x,mmax_var+0*mmax(1,:),'--k','LineWidth',2);
%         xlim([min(x)*.995 max(x)*1.0]);
        xlim([5 max(x)*1.0]);
        ylim([0 2.5])
        title('Government spending multiplier (peak)')
        legend([g2 g1],{'VAR','GMA'},'Location','NorthWest'), legend('boxoff')
        ylabel('Stimulus shock')
        set(gca,'ytick',0:.5:2.5)
    end
end


% MSUM:
mu_ur = 5.6;
x = mu_ur+INDICATOR;
figure,
for j=1:2
    error1_pos = msum(j,:) - msumlow(j,:);
    error2_pos = msumhigh(j,:) - msum(j,:);
    if j==1,
        subplot(212), hold on
        H= shadedErrorBar(x,msum(j,:),[error2_pos; error1_pos], 'r');
        g1 = plot(x,msum(j,:), 'r', 'LineWidth', 2);
        g2 = plot(x,msum_var+0*msum(1,:),'--k','LineWidth',2);
        xlim([min(x)*.995 max(x)*1.0]);
        legend([g2 g1],{'VAR','GMA'},'Location','Best'), legend('boxoff')     
        ylabel('Consolidation shock')
        ylim([0 2.5])
        xlabel('Unemployment rate')
        set(gca,'ytick',0:.5:2)
    end
    if j==2
        subplot(211), hold on
        H= shadedErrorBar(x,msum(j,:),[error2_pos; error1_pos], 'b')
        g1 = plot(x,msum(j,:), 'b', 'LineWidth', 2);
        g2 = plot(x,msum_var+0*msum(1,:),'--k','LineWidth',2);
        xlim([min(x)*.995 max(x)*1.0]);
        ylim([0 2.5])
        title('Government spending multiplier (integral)')
        legend([g2 g1],{'VAR','GMA'},'Location','Best'), legend('boxoff')
        ylabel('Stimulus shock')
        set(gca,'ytick',0:.5:2)
    end
end




