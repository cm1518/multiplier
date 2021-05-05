figure, scatter(squeeze(Msum(1,1,:)),squeeze(Msum(1,size(Msum,2),:)))
hold on, line([-.5 5],[-.5 5],'LineStyle','--','Color','r')
xlim([-.5 5])
ylim([-.5 5])

% Count # of points above 45degree line:
count=sum(squeeze(Msum(1,size(Msum,2),:))>squeeze(Msum(1,1,:)));

% Posterior proba that msum(UR high)>msum(UR low)
count/length(Msum)

    