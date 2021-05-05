function [ log_l_temp indicator_temp epsilon_temp ]=likelihood_increment(data, current_matrices,laggedmatrices,eps_initial, lag_length,indicator )
%computs the increment to the log likelihood function
LD=length(data);

%get reduced form errors
u=data;
for ll=1:lag_length
   u=u-squeeze(laggedmatrices(:,:,bin2dec(sprintf('%u',indicator(:,end-ll+1))),ll))*eps_initial(:,end-ll+1); 
  %use this for mean as well  
end

indicator_temp=zeros(LD,1);
for kk=1:LD
   indicator_temp(kk)=(u(kk)/currentmatrices
    
end
end

