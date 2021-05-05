
%input Ramey shocks to get their distribution over the business cycle
shocks=data(4,:);
r_BM=shocks;


% number of bins=N
N=14;

% close all
% indicator=setup.indicator+5.6;
indicator=setup.indicator;
indicator=indicator(lags+1:end);
[n,xout] = hist(indicator,N+1);

indicator=[indicator; nan(1,1)];

% % ATTENTION: need to forward indicator by 6 months (??)
% indicator0=indicator;
% indicator=[indicator(6:end) nan(1,6)];
%indicator=[nan(1,6) indicator(1:end-5) ];

% [row,col,v] = find((a(t,:)>r1(1,:)) & (a(t,:)<r2(1,:)));


bin_ind=N+1+zeros(length(indicator),1);
% 1) re-write indicator in bins
for j=1:length(indicator)
    k=1;
    while k<N+1
    if indicator(j)<xout(k)
        bin_ind(j)=k;
        k=N+1;
    else
        k=k+1;
    end
    end
end

    
% 2) assign the MP shocks to those bins
r_BM_bin=zeros(N+1,1);
r_BMpos_bin=zeros(N+1,1);
r_BMneg_bin=zeros(N+1,1);

r_BM_binv=zeros(N+1,1);
r_BMpos_binv=zeros(N+1,1);
r_BMneg_binv=zeros(N+1,1);


for j=1:length(r_BM)
    % I count the nb of MP shocks in each bin
    r_BM_bin((bin_ind(j)))=r_BM_bin(bin_ind(j))+(r_BM(j)~=0);
    r_BMpos_bin((bin_ind(j)))=r_BMpos_bin(bin_ind(j))+(r_BM(j)>0.0);
    r_BMneg_bin((bin_ind(j)))=r_BMneg_bin(bin_ind(j))+(r_BM(j)<0.0);
    
     % I sum the values of the MP shocks in each bin
    r_BM_binv((bin_ind(j)))=r_BM_binv(bin_ind(j))+sign(r_BM(j))*r_BM(j)*(r_BM(j)~=0);
    r_BMpos_binv((bin_ind(j)))=r_BMpos_binv(bin_ind(j))+r_BM(j)*(r_BM(j)>0.0);
    r_BMneg_binv((bin_ind(j)))=r_BMneg_binv(bin_ind(j))-r_BM(j)*(r_BM(j)<0.0);
end


%weighted:
%number of shocks:
[sigmap mup Ap]=mygaussfit(xout,r_BMpos_binv./sum(r_BMpos_binv) )
[sigman mun An]=mygaussfit(xout,r_BMneg_binv./sum(r_BMneg_binv) )

% Normalized-size shock distribution
figure,
plot(xout,[Ap*exp(-(xout-mup).^2/(2*sigmap^2)); ],'LineWidth',1.5)
hold on, plot(xout, An*exp(-(xout-mun).^2/(2*sigman^2)),'LineWidth',1.5,'Color',[0 .5 0])
% hold on, plot(12*xout,[r_BMpos_binv./sum(r_BMpos_binv) r_BMneg_binv./sum(r_BMneg_binv)],'--')
line([mean(setup.indicator) mean(setup.indicator)],[0 .14],'LineStyle',':','Color',[.4,.4,.4])
hold on, bar(xout,r_BMpos_binv./sum(r_BMpos_binv),1)
hold on, bar(xout,r_BMneg_binv./sum(r_BMneg_binv),1,'FaceColor',[0 .5 0])
legend('Expansionary shocks','Contractionary shocks','Location','Best'),legend('Boxoff')
xlabel('Unemployment rate')
% xlim([-3 3])
% ylim([0 .16])

% Actual-size shock distribution
figure,
bar(xout,r_BMpos_binv,1)
hold on, bar(xout,r_BMneg_binv,1,'FaceColor',[0 .5 0])
legend('Contractionary shocks','Expansionary shocks','Location','Best'),legend('Boxoff')
xlabel('Unemployment rate')

% stop
% 
% figure, plot(12*xout,[r_BM_bin./sum(r_BM_bin) r_BMpos_bin./sum(r_BMpos_bin) r_BMneg_bin./sum(r_BMneg_bin)])
% 
% %value of shocks:
% figure, plot(12*xout,[r_BM_binv./sum(r_BM_binv) r_BMpos_binv./sum(r_BMpos_binv) r_BMneg_binv./sum(r_BMneg_binv)])
% 
% %unweighted:
% figure, plot(12*xout,[r_BM_bin r_BMpos_bin r_BMneg_bin])
% figure, plot(12*xout,[r_BM_binv r_BMpos_binv r_BMneg_binv])
% 
% 
% [bandwidth,densityVAR_pos,xmeshVAR_pos]=kde(r_BMpos_bin,2^12);
% figure, plot(xmeshVAR_pos,densityVAR_pos,'r--','LineWidth',2)
% [bandwidth,densityVAR_neg,xmeshVAR_neg]=kde(r_BMneg_bin,2^12);
% hold on, plot(xmeshVAR_neg,densityVAR_neg,'b--','LineWidth',2)
% 
% 
% %alternative distribution (As in Tenreyo)
% %dist of shocks for espansins and recessions
% [pos bins]= hist(r_BM_pos,50);
% neg= hist(r_BM_neg,bins);
% figure, plot([pos]);
% hold on, plot(neg,'r');
% 
% 
% 
% break
% r_BM0=r_BM;
% % r_BM=tsmovavg(r_BM',1,12)';
% 
% r_BM_pos=nan(length(r_BM),1);
% r_BM_neg=r_BM_pos;
% for i=1:length(r_BM)
%     if r_BM(i)>=0
%         r_BM_pos(i)=r_BM(i);
%     end
%         if r_BM(i)<0
%         r_BM_neg(i)=r_BM(i);
%     end
% end
% 
% figure, scatter(indicator,r_BM_pos)
% hold on, scatter(indicator,r_BM_neg,'r')
% 
%  figure, plot([10*indicator' 10*r_BM])
%  
%  
%  figure, plot(tsmovavg(r_BM',1,12)')
%  
%  
 