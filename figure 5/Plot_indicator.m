close all
clear;
clc;


data_tp=xlsread(['RZ_Data_1890_2014.xlsx']);


setup.lags_shocks=45;

Tinit=1+ setup.lags_shocks;%setup.lags_shocks; %only start at 41 because setup.lags_shocks=40
    
data_tp=data_tp(Tinit-setup.lags_shocks:end,:);
yearn=data_tp(:,1);
nber=data_tp(:,6);

indicator_SUR=data_tp(:,50); %use Ramey's cyclical indicator (ur_cyc1)

% figure, plot(yearn, indicator_SUR)

figure1=figure,
hold on, plot(305+0*(-10:30),-10:30,'LineStyle',':','Color',[.8 .2 .2])
hold on, plot(1:length(indicator_SUR),2+0*indicator_SUR,'LineStyle','--','Color',[.6 .6 .6],'LineWidth',1)
hold on, plot(1:length(indicator_SUR),-1+0*indicator_SUR,'LineStyle','--','Color',[.6 .6 .6],'LineWidth',1)
plot(indicator_SUR,'LineWidth',2,'Color','k')
set(gca,'Xtick',1:40:length(indicator_SUR),'Xticklabel',num2str([yearn(1):10:yearn(end)]'),'FontSize',12)
set(gca,'Fontsize',8)
ylim([-5 20])

% Create textbox
annotation(figure1,'textbox',...
    [0.60 0.82 0.1 0.05],...
    'String',{'Post 1966'},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor',[1 0 0],...
    'Color',[.8 .2 .2]);

% Create arrow
annotation(figure1,'arrow',[0.62 0.67],...
    [0.816 0.816],'HeadLength',6,'HeadWidth',6,...
    'Color',[.8 .2 .2]);

% Create textbox
annotation(figure1,'textbox',[0.89 0.31 0.1 0.05],...
    'String',{'+2'},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor',[1 0 0],...
    'Color',[0.6 0.6 0.6]);

% Create textbox
annotation(figure1,'textbox',[0.90 0.205 0.1 0.05],...
    'String',{'-1'},...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor',[1 0 0],...
    'Color',[0.6 0.6 0.6]);




print('Figure5','-depsc') 
savefig('Figure5')
