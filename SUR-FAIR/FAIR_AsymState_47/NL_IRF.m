function [ Sigma_lagged ] = NL_IRF(a,b,c,beta,gamma,lags,ind)
%computes IRF pattern
lag_vec=0:lags-1;



Sigma_lagged=exp(-((lag_vec-b).^2)./c)...
    .*(a*(1+(gamma*exp(beta*ind)/(1+exp(beta*ind)))));



end

