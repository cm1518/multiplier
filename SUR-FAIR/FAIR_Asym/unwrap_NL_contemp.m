function [ Sigma] = unwrap_NL_contemp( params,epsilon_vec,setup )
%function returns contmeporaneous Sigma matrix - setup.size_obs by setup.size_obs 
%parametrs are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag


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
counter=counter+setup.size_obs^2;



Sigma=tril(SL_NL( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec,0,setup.size_obs ),-1) ;

for kk=1:setup.size_obs
   Sigma(kk,kk,1)=exp(alpha_diag(kk)*abs(epsilon_vec(kk,1))+beta_diag(kk)); 
    
end

end

