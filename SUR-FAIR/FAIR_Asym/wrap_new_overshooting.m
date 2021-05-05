function [ constant current_matrices lagged_matrices ] = wrap_BM( params, setup )


%first set of parameters are constants, then contemporaneous sigma_matrices
%order of parameters: first a then b, then c: all are vectors of size
%setup.size_obs^2 by 1. First , we have the parameters for negative shocks (first normal type term),
%then for positive (first normal type term) and then the same ordering for
%the second normal type term
%also imposes the identification restriction that that elements of |a_xxx| are larger
%than the correpsonding elements of |a_xxx2|. It is convenient to impose
%those restictions in this function even though it might be more natural to
%impose it when the propsal draw is made. 



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
ind_for_loop=ind_for_loop+setup.size_obs^2;


%parameters for negative shocks - second part
a_neg2=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
b_neg2=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
c_neg2=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;

%parameters for positive shocks - second part
a_pos2=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
b_pos2=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);
ind_for_loop=ind_for_loop+setup.size_obs^2;
c_pos2=repmat(reshape(params(ind_for_loop+1:ind_for_loop+setup.size_obs^2),setup.size_obs,setup.size_obs),1,setup.lags);



%imposing identification restriction that elements of |a_xxx| are larger
%than the correpsonding elements of |a_xxx2|

indneg=abs(a_neg)>abs(a_neg2);
indpos=abs(a_pos)>abs(a_pos2);

a_neg_new=indneg.*a_neg+(1-indneg).*a_neg2;
b_neg_new=indneg.*b_neg+(1-indneg).*b_neg2;
c_neg_new=indneg.*c_neg+(1-indneg).*c_neg2;

a_pos_new=indpos.*a_pos+(1-indpos).*a_pos2;
b_pos_new=indpos.*b_pos+(1-indpos).*b_pos2;
c_pos_new=indpos.*c_pos+(1-indpos).*c_pos2;


indneg=1-indneg;
indpos=1-indpos;

a_neg2_new=indneg.*a_neg+(1-indneg).*a_neg2;
b_neg2_new=indneg.*b_neg+(1-indneg).*b_neg2;
c_neg2_new=indneg.*c_neg+(1-indneg).*c_neg2;

a_pos2_new=indpos.*a_pos+(1-indpos).*a_pos2;
b_pos2_new=indpos.*b_pos+(1-indpos).*b_pos2;
c_pos2_new=indpos.*c_pos+(1-indpos).*c_pos2;

lagged_matrices(:,:,1,:)=reshape(a_neg_new.*exp((-(j_matrix-b_neg_new).^2)./c_neg_new)+a_neg2_new.*exp((-(j_matrix-b_neg2_new).^2)./c_neg2_new),setup.size_obs,setup.size_obs,setup.lags);
lagged_matrices(:,:,2,:)=reshape(a_pos_new.*exp((-(j_matrix-b_pos_new).^2)./c_pos_new)+a_pos2_new.*exp((-(j_matrix-b_pos2_new).^2)./c_pos2_new),setup.size_obs,setup.size_obs,setup.lags);

