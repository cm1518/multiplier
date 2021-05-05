function [ Sigma_lagged ] = SL_NL( alpha_gen,beta_gen,b_gen,c_gen,epsilon_vec,lag,size_obs, thresh_vec )
%computes one lagged Sigma matrix for the non-linear case

Sigma_lagged=exp(-((lag-reshape(b_gen,size_obs,size_obs)).^2)./reshape(c_gen,size_obs,size_obs))...
    .*(reshape(alpha_gen,size_obs,size_obs).*(abs(((epsilon_vec>=0).*(min(epsilon_vec,thresh_vec))+(epsilon_vec<0).*max(epsilon_vec,-thresh_vec))*ones(1,size_obs)))'+reshape(beta_gen,size_obs,size_obs));


end

