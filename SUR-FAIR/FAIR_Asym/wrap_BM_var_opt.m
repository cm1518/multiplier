function [ MA_matrices ] = wrap_BM_var_opt( params, setup )
%function that maps a parameter vector into the objects needed to evaluate
%the likelihood function

%transforming parameters back to the original parameter space (NOT
%NECESSARY FOR IRFS, FINAL DRAWS ALREADY IN UN-TRANSFORMED VARIABLES)
%[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );

params=[zeros(setup.size_obs,1);params];
constant=params(1:setup.size_obs);

current_matrices=zeros(setup.size_obs,setup.size_obs,2);

current_matrices(:,:,1)=ltvec(params(setup.size_obs+1:setup.size_obs+(setup.size_obs*(setup.size_obs+1))/2));
%current_matrices(:,:,2)=ltvec(params(setup.size_obs+(setup.size_obs*(setup.size_obs+1))/2+1:setup.size_obs+setup.size_obs*(setup.size_obs+1)/2+(setup.size_obs*(setup.size_obs+1))/2));
current_matrices(:,:,2)=current_matrices(:,:,1);
lagged_matrices=zeros(setup.size_obs,setup.size_obs,2,setup.lags);

MA_matrices=zeros(setup.size_obs,setup.size_obs,setup.lags+1);
MA_matrices(:,:,1)=current_matrices(:,:,1);
j_matrix=kron(1:1:setup.lags,ones(setup.size_obs,setup.size_obs));

ind_for_loop=setup.size_obs+setup.size_obs*(setup.size_obs+1)/2;

%parameters for negative shocks
a_neg=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
b_neg=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
c_neg=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;

%parameters for positive shocks
% a_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
% ind_for_loop=ind_for_loop+setup.size_obs^2;
% b_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
% ind_for_loop=ind_for_loop+setup.size_obs^2;
% c_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);

a_pos=a_neg;
b_pos=b_neg;
c_pos=c_neg;

lagged_matrices(:,:,1,:)=reshape(a_neg.*exp((-(j_matrix-b_neg).^2)./c_neg),setup.size_obs,setup.size_obs,setup.lags);
MA_matrices(:,:,2:end)=lagged_matrices(:,:,1,:);
lagged_matrices(:,:,2,:)=reshape(a_pos.*exp((-(j_matrix-b_pos).^2)./c_pos),setup.size_obs,setup.size_obs,setup.lags);

