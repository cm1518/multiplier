function [ log_l, indicator, epsilon, uvec] = likelihood(data,  params,setup)

%computes the log likelihood for the MA model
 eps_initial =setup.initial_eps; 
lag_length=setup.lags;


indicator=zeros(size(data,1),lag_length+size(data,2));
epsilon=zeros(size(data,1),lag_length+size(data,2));
%data is supposed to be a matrix with each column containing the
%observations at one point in time
uvec=zeros(size(data,1),size(data,2));

%careful about direction of epsilon/data!
%pass params as an argument to this function
%update epsilon! DONE

%setting options for optimization

options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',500, 'MaxFunEvals',5000);

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
alpha_gen_neg=reshape(alpha_gen_neg,setup.size_obs,setup.size_obs);
beta_gen_neg=reshape(beta_gen_neg,setup.size_obs,setup.size_obs);
b_gen_neg=reshape(b_gen_neg,setup.size_obs,setup.size_obs);
c_gen_neg=reshape(c_gen_neg,setup.size_obs,setup.size_obs);


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
alpha_gen_pos=reshape(alpha_gen_pos,setup.size_obs,setup.size_obs);
beta_gen_pos=reshape(beta_gen_pos,setup.size_obs,setup.size_obs);
b_gen_pos=reshape(b_gen_pos,setup.size_obs,setup.size_obs);
c_gen_pos=reshape(c_gen_pos,setup.size_obs,setup.size_obs);


epsilon(:,1:lag_length)=eps_initial;
%indicator equal to 0 if epsilon is negative or zero, 1 otherwise
log_l=0;
for tt=1:size(data,2)
%     if tt==1
%         epsilon(:,tt:tt+lag_length-1)
%     end

%solve for contemporaneous epsilons

eps_initial=[epsilon(:,tt:tt+lag_length-1) zeros(setup.size_obs,1)]; % we pad zeros on for the contmeporaneous epsilons

%first, get lagged sigma matrices
[ Sigma] = unwrap_NL_lagged2(  alpha_gen_neg,beta_gen_neg,b_gen_neg,c_gen_neg,alpha_gen_pos,beta_gen_pos,b_gen_pos,c_gen_pos,fliplr(eps_initial),setup );






condmean=zeros(setup.size_obs,1);

for ll=1:lag_length
   
   
   condmean=condmean+Sigma(:,:,ll)*eps_initial(:,end-ll); %note the reverse ordering of Sigma and epsilon
end
condmean=condmean+intercept;
u=data(:,tt)-condmean;
uvec(:,tt)=u;
temp=setup.store_responses(1:length(setup.index_restricted),1:length(setup.index_restricted),1);

epsilon_temp=zeros(size(data,1),1);
epsilon_temp(1:length(setup.index_restricted))=temp\u(1:length(setup.index_restricted));

temp2=setup.store_responses(end,1:length(setup.index_restricted),1)*epsilon_temp(1:length(setup.index_restricted));
if (((u(end)-temp2)/beta_diag_neg(end))<0)
    ind=0;
  epsilon_temp(end)=((u(end)-temp2)/beta_diag_neg(end));
else
    ind=1;
    epsilon_temp(end)=((u(end)-temp2)/beta_diag_pos(end));
end

% temp=setup.store_responses(:,:,1);
% for jj=length(setup.index_restricted)+1:setup.size_obs
%    temp(jj,jj)=beta_diag(jj); 
% end

%solve for restricted epsilons

%solve for unrestricted epsilons

% for kk=(length(setup.index_restricted)+1):setup.size_obs
%     [ sum_temp] = within_period_NL2( alpha_gen,beta_gen,b_gen,c_gen,epsilon_temp,setup,kk );
%     
% postopt2=@(shocks) NLcontemp_abs(shocks,alpha_diag(kk),beta_diag(kk),u(kk)-sum_temp);
%   
% 
% % if tt==1
% %     start=0;
% % else
% %    start=epsilon_temp(kk-1); 
% % end
% 
% [epsilon_temp(kk),~]=fzero(postopt2,0);
% %[epsilon_temp(kk),~]=fsolve(postopt2,0,options);
% 
% 
% 
% end

% [ sigma_temp] = unwrap_NL_contemp( params,epsilon_vec,setup );
% 
% covariance=sigma_temp*sigma_temp';

Sigma_contemp=setup.store_responses(:,:,1);
Sigma_contemp(end,end)=(1-ind)*beta_diag_neg(end)+ind*beta_diag_pos(end);
covariance=Sigma_contemp*Sigma_contemp';
epsilon(:,tt+lag_length)=epsilon_temp;
%log_l_temp_vec=zeros(setup.size_obs,1);
%for jj=1:setup.size_obs
log_l_temp=-(setup.size_obs/2)*log(2*pi)-.5*log(det(covariance))-.5*((Sigma_contemp*epsilon_temp)'/covariance)*(Sigma_contemp*epsilon_temp);
%end



log_l=log_l+log_l_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
if (sum(epsilon(3,:)==0)-setup.lags)>2
    log_l=-inf;
end