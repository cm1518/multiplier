function [ sum] = within_period_NL( params,epsilon_vec,setup,order_of_var )
%function returns sum of elements below main diagonal for contemporanous
%matrix

if order_of_var==1
    
    sum=0;
    
elseif order_of_var>1
%unwrapping parameter vector
intercept=params(1:setup.size_obs);
counter=setup.size_obs+1;

alpha_diag=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

beta_diag=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

alpha_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

b_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

c_gen=params(counter:counter+setup.size_obs^2-1);





Sigma=tril(SL_NL( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec,0,setup.size_obs ),-1) ;
 


sum_temp=Sigma*epsilon_vec;
sum=sum_temp(order_of_var);
end
end
