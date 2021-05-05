function [ sum] = within_period_NL2( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec,setup,order_of_var )
%function returns sum of elements below main diagonal for contemporanous
%matrix

if order_of_var==1
    
    sum=0;
    
elseif order_of_var>1
%unwrapping parameter vector




Sigma=tril(SL_NL_2( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec,0,setup.size_obs ),-1) ;
 
Sigma(:,setup.index_restricted)=setup.store_responses(:,setup.index_restricted,1);
Sigma=tril(Sigma,-1);


sum_temp=Sigma*epsilon_vec;
sum=sum_temp(order_of_var);
end
end
