function [ xestimate1, functionvalue1] = first_pass_estimation( setup )
%estimates the coefficients on the indicator keeping all other parameters
%fixed


data=load(setup.data);
data=data.data;
warning off;
%options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-8,'MaxIter',1e5, 'MaxFunEvals',1e6);
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxIter',500, 'MaxFunEvals',50000);


 postopt2=@(params) (-likelihood_wrap_2( params,setup,data ));
 
x0=0*ones(setup.size_obs*2,1); 
 
 [xestimate1,functionvalue1]=fminsearch(postopt2,x0,options);
%[xestimate1,functionvalue1]=fminunc(postopt2,x0,options);
