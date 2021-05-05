function [ Sigma_lagged ] = SL_NL_3(alpha_gen_neg, beta_gen_neg,b_gen_neg,c_gen_neg,alpha_gen_pos,beta_gen_pos,b_gen_pos,c_gen_pos,epsilon_vec,lag,size_obs,indicator,mean_ind )
%computes one lagged Sigma matrix for the non-linear case
ind=epsilon_vec<0*ones(size(alpha_gen_neg));


	
alpha_gen=alpha_gen_neg.*(ind)+alpha_gen_pos.*((1-ind));
beta_gen=repmat(beta_gen_neg,1,size(alpha_gen,2)).*(ind)+repmat(beta_gen_pos,1,size(alpha_gen,2)).*((1-ind));



b_gen=b_gen_neg.*(ind)+b_gen_pos.*((1-ind));
c_gen=c_gen_neg.*(ind)+c_gen_pos.*((1-ind));

Sigma_lagged=exp(-((lag-b_gen).^2)./c_gen)...
    .*(alpha_gen.*(1+beta_gen*(indicator-mean_ind)));



	
end

