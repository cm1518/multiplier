function [ Sigma, intercept] = unwrap_NL_IRF_peak( params,epsilon_vec,setup,indicator_vec )
%function returns intercept and array of Sigma matrices - setup.size_obs by setup.size_obs by setup.lag_length 
%with Sigma_0 ordered first
%parametrs are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag
[ params ] = params_mod( params,setup );

%unwrapping parameters

intercept=params(1:setup.size_obs);
counter=setup.size_obs+1;

beta_diag_neg=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

alpha_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

b_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

c_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;
gamma_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;
alpha_gen_neg=reshape(alpha_gen_neg,setup.size_obs,setup.size_obs);
beta_gen_neg=reshape(beta_gen_neg,setup.size_obs,setup.size_obs);
b_gen_neg=reshape(b_gen_neg,setup.size_obs,setup.size_obs);
c_gen_neg=reshape(c_gen_neg,setup.size_obs,setup.size_obs);
gamma_gen_neg=reshape(gamma_gen_neg,setup.size_obs,setup.size_obs);

beta_diag_pos=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

alpha_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

b_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

c_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;
gamma_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;
alpha_gen_pos=reshape(alpha_gen_pos,setup.size_obs,setup.size_obs);
beta_gen_pos=reshape(beta_gen_pos,setup.size_obs,setup.size_obs);
b_gen_pos=reshape(b_gen_pos,setup.size_obs,setup.size_obs);
c_gen_pos=reshape(c_gen_pos,setup.size_obs,setup.size_obs);
gamma_gen_pos=reshape(gamma_gen_pos,setup.size_obs,setup.size_obs);


Sigma=zeros(setup.size_obs,setup.size_obs,setup.lags+1);


Sigma(:,:,1)=setup.store_responses(:,:,1) ;

for jj=length(setup.index_restricted)+1:setup.size_obs
   Sigma(jj,jj,1)=beta_diag_neg(jj)*(epsilon_vec(jj,1)<0)+beta_diag_pos(jj)*(epsilon_vec(jj,1)>=0); 
end

for jj=1:setup.lags

    Sigma(:,:,jj+1)= SL_NL_2_exp( alpha_gen_neg,beta_gen_neg,b_gen_neg,c_gen_neg,gamma_gen_neg,alpha_gen_pos,beta_gen_pos,b_gen_pos,c_gen_pos,gamma_gen_pos,epsilon_vec(:,jj+1),jj,setup.size_obs,setup.threshold_vec,indicator_vec(end-jj+1));

end
Sigma(:,setup.index_restricted,:)=setup.store_responses(:,setup.index_restricted,1:end);

end

