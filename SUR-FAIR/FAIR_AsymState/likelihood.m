function [ log_l,uvec] = likelihood(data,  params,setup)

%computes the log likelihood for the SUR model

%RHS data common to all dependent variables
data_LHS=data(1:setup.dep_variables,:);
data_RHS_linear=data(setup.dep_variables+1:setup.dep_variables+1+setup.lags_endo+setup.lags_exo*setup.dim_exo,:); %RHS variables except for shocks
data_RHS_shocks=data(setup.dep_variables+1+setup.lags_endo+setup.lags_exo*setup.dim_exo+1:end,:); %shocks

%vector of errors
uvec=zeros(setup.dep_variables,size(data,2));
 counter=1;
 for jj=1:setup.dep_variables

%unwrapping parameters

intercept{jj}=params(counter);
counter=1+counter;



lagged_y_param{jj}=params(counter:counter+setup.lags_endo-1);
counter=counter+setup.lags_endo;

lagged_exo_param{jj}=params(counter:counter+setup.lags_exo*setup.dim_exo-1); %I assume that parameters here are ordered by variable first, then by lag
counter=counter+setup.lags_exo*setup.dim_exo;


if jj==1
counter_linear=counter; %number of 'linear' (except AR coefficients) plus 1, same for all dep. variables
end

NL_a_neg{jj}=zeros(setup.num_Gaussian{jj},1);

NL_b_neg{jj}=zeros(setup.num_Gaussian{jj},1);

NL_c_neg{jj}=zeros(setup.num_Gaussian{jj},1);



NL_a_pos{jj}=zeros(setup.num_Gaussian{jj},1);

NL_b_pos{jj}=zeros(setup.num_Gaussian{jj},1);

NL_c_pos{jj}=zeros(setup.num_Gaussian{jj},1);

for ll=1:setup.num_Gaussian{jj}
  NL_a_neg{jj}(ll)=params(counter);
  counter=counter+1;
  NL_b_neg{jj}(ll)=params(counter);
  counter=counter+1;
  NL_c_neg{jj}(ll)=params(counter);
  counter=counter+1;
  NL_a_pos{jj}(ll)=params(counter);
  counter=counter+1;
  NL_b_pos{jj}(ll)=params(counter);
  counter=counter+1;
  NL_c_pos{jj}(ll)=params(counter);
  counter=counter+1;
end


gamma_neg{jj}=params(counter);
  counter=counter+1;
  
 
  

beta_neg{jj}=params(counter);
  counter=counter+1;
  
  
   gamma_pos{jj}=params(counter);
  counter=counter+1;
beta_pos{jj}=params(counter);
  counter=counter+1;



  





end
sigma_u=ltvec(params(counter:counter+setup.dep_variables*(setup.dep_variables+1)/2-1));
sigma_u=sigma_u+sigma_u'; %this makes sure the off diagonal elements are correct, but the diagonal elements are now twice the
%variances. Correct this in the next step.

%first, save variance terms
standard_devs=sqrt(.5*diag(sigma_u));

%then create correlation matrix
sigma_u(eye(size(sigma_u))~=0)=1;

%finally compute covariance matrix

sigma_u=diag(standard_devs)*sigma_u*diag(standard_devs);

ar_terms=reshape(params(counter+setup.dep_variables*(setup.dep_variables+1)/2:end),setup.dep_variables,setup.dep_variables);



lagged_u=setup.ar_initial;
log_l=0;
cond_mean_rec=zeros(setup.dep_variables,1);
for tt=1:size(data,2)
%note that uvec(tt) is now the iid error in the AR modell for the residual,
%whereas lagged_u are the lagged autocorrelated residuals

ind=fliplr(setup.indicator(tt:tt+setup.lags_shocks-1))';
for jj=1:setup.dep_variables

Sigma_lagged_neg{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
for gg=1:setup.num_Gaussian{jj}
    [ Sigma_lagged_neg{jj}(gg,:) ] = NL_IRF(NL_a_neg{jj}(gg),NL_b_neg{jj}(gg),NL_c_neg{jj}(gg),gamma_neg{jj},beta_neg{jj},setup.lags_shocks,ind); 
    
end

Sigma_lagged_pos{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
for gg=1:setup.num_Gaussian{jj}
    [ Sigma_lagged_pos{jj}(gg,:) ] = NL_IRF(NL_a_pos{jj}(gg),NL_b_pos{jj}(gg),NL_c_pos{jj}(gg),gamma_pos{jj},beta_pos{jj},setup.lags_shocks,ind); 
    
end



cond_mean_rec(jj)=[intercept{jj} lagged_y_param{jj}' lagged_exo_param{jj}']*data_RHS_linear(:,tt)+ar_terms(jj,:)*lagged_u;

cond_mean_rec(jj)=cond_mean_rec(jj)+sum(Sigma_lagged_neg{jj}*((data_RHS_shocks(:,tt)<0).*data_RHS_shocks(:,tt)))+sum(Sigma_lagged_pos{jj}*((data_RHS_shocks(:,tt)>=0).*data_RHS_shocks(:,tt)));

end
%calcualte iid error term
uvec(:,tt)=data_LHS(:,tt)-cond_mean_rec;

%calculate lagged AR error terms


lagged_u=[uvec(:,tt)+ar_terms*lagged_u];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%log_l_temp_vec=zeros(setup.size_obs,1);
%for jj=1:setup.size_obs
log_l_temp=-(setup.dep_variables/2)*log(2*pi)-.5*logdet(sigma_u)-.5*(uvec(:,tt)'/sigma_u*uvec(:,tt));
%end



log_l=log_l+log_l_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
