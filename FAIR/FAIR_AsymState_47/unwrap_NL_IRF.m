function [ Sigma, intercept] = unwrap_NL_IRF( params,epsilon_vec,setup,indicator_vec,mean_ind)
%function returns intercept and array of Sigma matrices - setup.size_obs by setup.size_obs by setup.lag_length 
%with Sigma_0 ordered first
%parameters are ordered as follows: [intercepts;alpha_diags;beta_diags;alpha_gen;beta_gen;b_gen;c_gen]
%epsilons are ordered from most recent to epsilon with largest lag
[ params ] = params_mod( params,setup );

%unwrapping parameters


Sigma=zeros(setup.size_obs,setup.size_obs,setup.lags+1);


Sigma(:,:,1)=setup.store_responses(:,:,1) ;



if epsilon_vec(setup.index_unrestricted,1)<0
temp_neg=setup.store_responses(:,:,1);
temp2=sum(params.alpha_gen_neg.*repmat((1+(params.beta_gen_neg*(indicator_vec(end)-mean_ind))),1,size(params.alpha_gen_neg,2)),2);
temp_neg(setup.index_unrestricted:end,setup.index_unrestricted)=params.r_neg.*[temp2(setup.index_unrestricted:end,:)];

Sigma(:,:,1)=temp_neg;
elseif epsilon_vec(setup.index_unrestricted,1)>=0
temp_pos=setup.store_responses(:,:,1);
temp2=sum(params.alpha_gen_pos.*repmat((1+(params.beta_gen_pos*(indicator_vec(end)-mean_ind))),1,size(params.alpha_gen_pos,2)),2);

temp_pos(setup.index_unrestricted:end,setup.index_unrestricted)=params.r_pos.*[temp2(setup.index_unrestricted:end,:)];
Sigma(:,:,1)=temp_pos;

end




 



%changed the indicator index here - latest entry must be dated end-1 because the last indicator is used to compute contemporaneous response

for jj=1:setup.lags
   % for nn=1:max(setup.num_gaussian)
  % Sigma(:,setup.index_unrestricted,jj+1)=SL_NL_2( params.alpha_gen_neg(:,nn),params.beta_gen_neg(:,nn),params.b_gen_neg(:,nn),params.c_gen_neg(:,nn),params.alpha_gen_pos(:,nn),params.beta_gen_pos(:,nn),params.b_gen_pos(:,nn),params.c_gen_pos(:,nn),epsilon_vec(setup.index_unrestricted,jj+1),jj,setup.size_obs,indicator_vec(end-jj),mean_ind) +Sigma(:,setup.index_unrestricted,jj+1);
   % end
   jj;
   Sigma(:,setup.index_unrestricted,jj+1)=sum(SL_NL_3( params.alpha_gen_neg(:,:),params.beta_gen_neg(:,:),params.b_gen_neg(:,:),params.c_gen_neg(:,:),params.alpha_gen_pos(:,:),params.beta_gen_pos(:,:),params.b_gen_pos(:,:),params.c_gen_pos(:,:),epsilon_vec(setup.index_unrestricted,jj+1),jj,setup.size_obs,indicator_vec(end-jj),mean_ind),2);
  
end


Sigma(:,setup.index_restricted,:)=setup.store_responses(:,setup.index_restricted,1:end);
intercept=params.intercepts;
end

