clear;
close all;
clc;
%Plot distribution of shocks:

%Use AG shocks:
load gorod_ext.mat
shocksAG=Gshock_Ghat; %Shocks using professional forecasts as control

%USE 1890-2014 DATA from Ramey
data=xlsread(['RZ_Data.xlsx']);


Yn=data(:,3);
news=data(:,8);
%re-scale news as share of GDP
news(2:end)=news(2:end)./Yn(1:end-1);
ramey=100*news;
% load shock series
shocksRamey=ramey; %Ramey news shock


[i]=find(((shocksRamey)==0));

shocksRamey=shocksRamey;
shocksRamey(i)=NaN;

[i]=find(((shocksAG)==0));
shocksAG(i)=NaN;

figure, subplot(211), (hist(shocksAG/nanstd(shocksAG),10))
ylim([0 65]), xlim([-4 6]), title('Shocks distribution: Recursive identification')

subplot(212), (hist(shocksRamey/nanstd(shocksRamey),15))
ylim([0 65]), xlim([-4 6]), title('Shocks distribution: Narrative identification')


print('Figure8','-depsc') 
savefig('Figure8')


