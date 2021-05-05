function [ Sigma, intercept] = unwrap_NL_IRF( params,epsilon_vec,setup,indicator_vec)
%function returns intercept and array of Sigma matrices - setup.size_obs by setup.size_obs by setup.lag_length 
%with Sigma_0 ordered first
%parametrs are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag
[ params ] = params_mod( params,setup );

%unwrapping parameters


Sigma=zeros(setup.size_obs,setup.size_obs,setup.lags+1);


Sigma(:,:,1)=setup.store_responses(:,:,1) ;



Sigma(setup.index_unrestricted:end,setup.index_unrestricted,1)=params.beta_diag_neg;





for jj=1:setup.lags
    for nn=1:max(setup.num_gaussian)
   Sigma(:,setup.index_unrestricted,jj+1)=SL_NL_2( params.beta_gen_neg(:,nn),params.b_gen_neg(:,nn),params.c_gen_neg(:,nn),jj) +Sigma(:,setup.index_unrestricted,jj+1);
    end
end


Sigma(:,setup.index_restricted,:)=setup.store_responses(:,setup.index_restricted,1:end);
intercept=params.intercepts;
end

