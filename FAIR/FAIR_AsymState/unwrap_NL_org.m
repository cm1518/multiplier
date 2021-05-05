function [ Sigma, intercept] = unwrap_NL_org( params,epsilon_vec,setup )
%function returns intercept and array of Sigma matrices - setup.size_obs by setup.size_obs by setup.lag_length 
%with Sigma_0 ordered first
%parametrs are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag


%unwrapping parameter vector
intercept=params(1:setup.size_obs);
counter=setup.size_obs+1;

beta_diag=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

alpha_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

b_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

c_gen=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

Sigma=zeros(setup.size_obs,setup.size_obs,setup.lags+1);


Sigma(:,:,1)=setup.store_responses(:,:,1) ;

for jj=length(setup.index_restricted)+1:setup.size_obs
   Sigma(jj,jj,1)=beta_diag(jj); 
end

for jj=1:setup.lags
   Sigma(:,:,jj+1)=SL_NL( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec(:,jj+1),jj,setup.size_obs,setup.threshold_vec ) ;
    
end
Sigma(:,setup.index_restricted,:)=setup.store_responses(:,setup.index_restricted,1:end);

end

