T=100;

rho=0.95;

data=zeros(1,T+1);

for tt=1:T
   data(tt+1)=rho*data(tt)+randn; 
    
end

data=data(2:end);

save('data_file.mat','data');