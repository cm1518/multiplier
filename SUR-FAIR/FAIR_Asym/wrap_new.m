function [ constant current_matrices lagged_matrices ] = wrap_BM( params, setup )


%first set of parameters are constants, then contemporaneous sigma_matrices
%order of parameters: first a then b, then c: all are vectors of size
%setup.size_obs^2 by 1. First , we have the parameters for negative shocks,
%then for positive



[ params ] = inv_transform( params,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );


constant=params(1:setup.size_obs);

current_matrices=zeros(setup.size_obs,setup.size_obs,2);

current_matrices(:,:,1)=ltvec(params(setup.size_obs+1:setup.size_obs+(setup.size_obs*(setup.size_obs+1))/2));
current_matrices(:,:,2)=ltvec(params(setup.size_obs+(setup.size_obs*(setup.size_obs+1))/2+1:setup.size_obs+setup.size_obs*(setup.size_obs+1)/2+(setup.size_obs*(setup.size_obs+1))/2));

lagged_matrices=zeros(setup.size_obs,setup.size_obs,2,setup.lags);


j_matrix=kron(1:1:setup.lags,ones(setup.size_obs,setup.size_obs));

ind_for_loop=setup.size_obs+setup.size_obs*(setup.size_obs+1)/2+(setup.size_obs*(setup.size_obs+1))/2;

%parameters for negative shocks
a_neg=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
b_neg=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
c_neg=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;

%parameters for positive shocks
a_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
b_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
c_pos=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);



lagged_matrices(:,:,1,:)=a_neg.*exp((-(j_matrix-b_neg).^2)./c_neg);
lagged_matrices(:,:,2,:)=a_pos.*exp((-(j_matrix-b_pos).^2)./c_pos);

