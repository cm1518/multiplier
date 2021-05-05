function [ log_l_temp indicator_temp epsilon_temp u]=likelihood_increment(data,intercept, current_matrices,laggedmatrices,eps_initial,  lag_length,indicator )
%computs the increment to the log likelihood function
LD=length(data);

%get reduced form errors
condmean=zeros(LD,1);

for ll=1:lag_length
   laggedtemp=zeros(LD,LD);
    for mm=1:LD
       laggedtemp(:,mm)=(laggedmatrices(:,mm,indicator(mm,end-ll+1)+1,ll));
    end
   
   condmean=condmean+laggedtemp*eps_initial(:,end-ll+1);
end
condmean=condmean+intercept;
u=data-condmean;



epsilon_temp=zeros(LD,1);
indicator_temp=zeros(LD,1);
sigma_temp=zeros(LD,LD);
for kk=1:LD
    temp=(u(kk)-sigma_temp(kk,:)*epsilon_temp)/current_matrices(kk,kk,2);
   indicator_temp(kk)=(temp)>0;
   if temp>0
       epsilon_temp(kk)=temp;
       sigma_temp(:,kk)=current_matrices(:,kk,2);
   else
      
      
       epsilon_temp(kk)=((u(kk)-sigma_temp(kk,:)*epsilon_temp)/current_matrices(kk,kk,1));
      sigma_temp(:,kk)=current_matrices(:,kk,1); 
   end
    
end


%log_l_temp2=log(mvnpdf(data,condmean,sigma_temp*sigma_temp'));
%CHANGE LOG DET CALCULATION!!!!
covariance=sigma_temp*sigma_temp';

% Using the residuals to build the likelihood:
log_l_temp=-(LD/2)*log(2*pi)-.5*log(det(covariance))-.5*(((data-condmean)')/covariance)*(data-condmean);

% Alternatively, using the structural shocks to build the likelihood
% log_l_temp_vec=zeros(3,1);
% for jj=1:3
% log_l_temp_vec(jj)=-(1/2)*log(2*pi)-.5*((epsilon_temp(jj))')*(epsilon_temp(jj));
% end
% log_l_temp=sum(log_l_temp_vec);


%log_l_temp-log_l_temp2
end

