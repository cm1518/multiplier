T=100;

rho=0.1;
rho_2=0.6;

data=zeros(2,T+1);

for tt=1:T
   data(1,tt+1)=rho*data(1,tt)+randn; 
    data(2,tt+1)=rho_2*data(2,tt)+randn;
end

data=data(:,2:end);

save('data_file.mat','data');