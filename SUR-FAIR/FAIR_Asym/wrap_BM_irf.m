function [ constant current_matrices lagged_matrices ] = wrap_BM_irf( params, setup )
%function that maps a parameter vector into the objects needed to evaluate
%the likelihood function

%transforming parameters back to the original parameter space (NOT
%NECESSARY FOR IRFS, FINAL DRAWS ALREADY IN UN-TRANSFORMED VARIABLES)
%[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );



constant=params(1:setup.size_obs);

current_matrices=zeros(setup.size_obs,setup.size_obs,2);
ind_for_loop=setup.size_obs;
current_matrices(:,:,1)=setup.store_responses(:,:,1);
for jj=1:length(setup.index_asymmetric_vars)
	hh=setup.index_asymmetric_vars(jj);
current_matrices(hh:end,hh,1)=reshape(params(ind_for_loop+1:ind_for_loop+1+setup.size_obs-hh),setup.size_obs-hh+1,1);
ind_for_loop=ind_for_loop+1+setup.size_obs-hh;
end


current_matrices(:,:,2)=current_matrices(:,:,1);
for jj=1:length(setup.index_asymmetric_vars)
	hh=setup.index_asymmetric_vars(jj);
current_matrices(hh:end,hh,2)=reshape(params(ind_for_loop+1:ind_for_loop+1+setup.size_obs-hh),setup.size_obs-hh+1,1);
ind_for_loop=ind_for_loop+1+setup.size_obs-hh;
end


lagged_matrices=zeros(setup.size_obs,setup.size_obs,2,setup.lags);


j_matrix=kron(1:1:setup.lags,ones(setup.size_obs,setup.size_obs));


%parameters for negative shocks
a_mat_neg=setup.a_VAR;
a_mat_neg(:,setup.index_asymmetric_vars)=(reshape(params(ind_for_loop+1:ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs),setup.size_obs,length(setup.index_asymmetric_vars)));
a_neg=repmat(a_mat_neg,1,setup.lags);
ind_for_loop=ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs;

b_mat_neg=setup.b_VAR;
b_mat_neg(:,setup.index_asymmetric_vars)=(reshape(params(ind_for_loop+1:ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs),setup.size_obs,length(setup.index_asymmetric_vars)));
b_neg=repmat(b_mat_neg,1,setup.lags);
ind_for_loop=ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs;

c_mat_neg=setup.c_VAR;
c_mat_neg(:,setup.index_asymmetric_vars)=(reshape(params(ind_for_loop+1:ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs),setup.size_obs,length(setup.index_asymmetric_vars)));
c_neg=repmat(c_mat_neg,1,setup.lags);
ind_for_loop=ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs;

%parameters for positive shocks
a_mat_pos=a_mat_neg;
a_mat_pos(:,setup.index_asymmetric_vars)=(reshape(params(ind_for_loop+1:ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs),setup.size_obs,length(setup.index_asymmetric_vars)));
a_pos=repmat(a_mat_pos,1,setup.lags);
ind_for_loop=ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs;


b_mat_pos=b_mat_neg;
b_mat_pos(:,setup.index_asymmetric_vars)=(reshape(params(ind_for_loop+1:ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs),setup.size_obs,length(setup.index_asymmetric_vars)));
b_pos=repmat(b_mat_pos,1,setup.lags);
ind_for_loop=ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs;

c_mat_pos=c_mat_neg;
c_mat_pos(:,setup.index_asymmetric_vars)=(reshape(params(ind_for_loop+1:ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs),setup.size_obs,length(setup.index_asymmetric_vars)));
c_pos=repmat(c_mat_pos,1,setup.lags);
ind_for_loop=ind_for_loop+length(setup.index_asymmetric_vars)*setup.size_obs;





lagged_matrices(:,:,1,:)=reshape(a_neg.*exp((-(j_matrix-b_neg).^2)./c_neg),setup.size_obs,setup.size_obs,setup.lags);
lagged_matrices(:,:,2,:)=reshape(a_pos.*exp((-(j_matrix-b_pos).^2)./c_pos),setup.size_obs,setup.size_obs,setup.lags);


lagged_matrices(:,setup.index_symmetric_vars,1,:)=setup.store_responses(:,setup.index_symmetric_vars,2:end);
lagged_matrices(:,setup.index_symmetric_vars,2,:)=setup.store_responses(:,setup.index_symmetric_vars,2:end);

 %lagged_matrices(:,:,1,:)=reshape(c_neg,setup.size_obs,setup.size_obs,setup.lags);
 %lagged_matrices(:,:,2,:)=reshape(c_pos,setup.size_obs,setup.size_obs,setup.lags);
