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

%options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxIter',500, 'MaxFunEvals',5000);

%build X matrix of RHS trend variables
X=[];
for jj=1:setup.polynomials+1
   % X=[X [0:setup.sample_size-1]'.^(jj-1)];
    X=[X (2*[0:setup.sample_size-1]'/(setup.sample_size-1)-1).^(jj-1)]; %rescale polynomials to make algorithm more stable
    % see https://www.rssd.esa.int/SP/LISAPATHFINDER/docs/Data_Analysis/DA_Nine/Annex7b.pdf
end


[ params ] = params_mod( params,setup );
%unwrapping parameters

mean_ind=mean(setup.indicator);
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
[ Sigma] = unwrap_NL_lagged2(  params.beta_gen_neg,params.b_gen_neg,params.c_gen_neg,fliplr(eps_initial),setup);






condmean=zeros(setup.size_obs,1);

for ll=1:lag_length
   
   
   condmean=condmean+Sigma(:,:,ll)*eps_initial(:,end-ll); %note the reverse ordering of Sigma and epsilon
end
condmean=condmean+[params.intercepts params.poly_coefficients]*X(tt,:)';
u=data(:,tt)-condmean;
uvec(:,tt)=u;
temp_neg=setup.store_responses(:,:,1);
temp_neg(setup.index_unrestricted:end,setup.index_unrestricted)=params.beta_diag_neg;



epsilon_temp=temp_neg\u;

    ind=0;
   temp_contemp=temp_neg;



Sigma_contemp=temp_contemp;
covariance=Sigma_contemp*Sigma_contemp';
epsilon(:,tt+lag_length)=epsilon_temp;
%log_l_temp_vec=zeros(setup.size_obs,1);
%for jj=1:setup.size_obs
log_l_temp=-(setup.size_obs/2)*log(2*pi)-.5*log(det(covariance))-.5*((Sigma_contemp*epsilon_temp)'/covariance)*(Sigma_contemp*epsilon_temp);
%end



log_l=log_l+log_l_temp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
