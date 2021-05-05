%set draw to whatever draw you want to use: median, mean or use a subset of
%draws and loop over the calculations below to get error bands

ind=1;
for size_of_shock=-2:.1:2;
    ind
IRFs=zeros(setup.size_obs,setup.lags+1);

for jj=1:setup.lags+1
   epsilon=[zeros(setup.size_obs,jj-1) [zeros(setup.size_obs-1,1);size_of_shock] zeros(setup.size_obs,setup.lags+1-jj)]; 
    [ Sigma, intercept] = unwrap_NL_IRF( draw,epsilon,setup );
    IRFs(:,jj)=Sigma(:,end,jj);
end
peak(:,ind)=max(abs(IRFs),[],2);
% figure;
% 
% for jj=1:setup.size_obs
% subplot(2,2,jj)
% plot(1:setup.lags+1,IRFs(jj,:),'Linewidth',2)
% title(sprintf('scaling, factor,shock size = %i', size_of_shock))
% end
% print -depsc


% for jj=1:setup.lags+1
%    epsilon=[zeros(setup.size_obs,jj-1) [zeros(setup.size_obs-1,1);size_of_shock] zeros(setup.size_obs,setup.lags+1-jj)]; 
%     [ Sigma, intercept] = unwrap_NL_IRF( draw,epsilon,setup );
%     IRFs(:,jj)=Sigma(:,end,jj)*size_of_shock;
% end

% figure;
% 
% for jj=1:setup.size_obs
% subplot(2,2,jj)
% plot(1:setup.lags+1,IRFs(jj,:),'Linewidth',2)
% title(sprintf('IRF, shock size = %i', size_of_shock))
% end
% 
% print -depsc
ind=ind+1;

end
figure
for jj=1:3
    subplot(2,2,jj)
    plot(-2:.1:2,peak(jj,:))
    title('peak of \Sigma')
end
print -depsc



ind=1;
for size_of_shock=-2:.1:2;
    ind
IRFs=zeros(setup.size_obs,setup.lags+1);

for jj=1:setup.lags+1
   epsilon=[zeros(setup.size_obs,jj-1) [zeros(setup.size_obs-1,1);size_of_shock] zeros(setup.size_obs,setup.lags+1-jj)]; 
    [ Sigma, intercept] = unwrap_NL_IRF( draw,epsilon,setup );
    IRFs(:,jj)=Sigma(:,end,jj)*size_of_shock;
end
peak(:,ind)=max(abs(IRFs),[],2);
% figure;
% 
% for jj=1:setup.size_obs
% subplot(2,2,jj)
% plot(1:setup.lags+1,IRFs(jj,:),'Linewidth',2)
% title(sprintf('scaling, factor,shock size = %i', size_of_shock))
% end
% print -depsc


% for jj=1:setup.lags+1
%    epsilon=[zeros(setup.size_obs,jj-1) [zeros(setup.size_obs-1,1);size_of_shock] zeros(setup.size_obs,setup.lags+1-jj)]; 
%     [ Sigma, intercept] = unwrap_NL_IRF( draw,epsilon,setup );
%     IRFs(:,jj)=Sigma(:,end,jj)*size_of_shock;
% end

% figure;
% 
% for jj=1:setup.size_obs
% subplot(2,2,jj)
% plot(1:setup.lags+1,IRFs(jj,:),'Linewidth',2)
% title(sprintf('IRF, shock size = %i', size_of_shock))
% end
% 
% print -depsc
ind=ind+1;

end
figure
for jj=1:3
    subplot(2,2,jj)
    plot(-2:.1:2,peak(jj,:))
    title('peak of IRF')
end
print -depsc