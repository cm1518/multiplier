function [ Sigma_lagged ] = SL_NL_2(alpha_gen_neg, beta_gen_neg,b_gen_neg,c_gen_neg,alpha_gen_pos,beta_gen_pos,b_gen_pos,c_gen_pos,epsilon_vec,lag,size_obs,indicator,mean_ind )
%computes one lagged Sigma matrix for the non-linear case
ind=epsilon_vec<0;




beta_gen=beta_gen_neg.*(ind*ones(1,size_obs))'+beta_gen_pos.*((1-ind)*ones(1,size_obs))';



b_gen=b_gen_neg.*(ind*ones(1,size_obs))'+b_gen_pos.*((1-ind)*ones(1,size_obs))';
c_gen=c_gen_neg.*(ind*ones(1,size_obs))'+c_gen_pos.*((1-ind)*ones(1,size_obs))';
alpha_gen=alpha_gen_neg.*(ind*ones(1,size_obs))'+alpha_gen_pos.*((1-ind)*ones(1,size_obs))';


Sigma_lagged=exp(-((lag-b_gen).^2)./c_gen)...
    .*(alpha_gen.*(1+beta_gen*(indicator-mean_ind)));



end

