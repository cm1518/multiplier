function [ Sigma_lagged ] = NL_IRF(a,b,c,lags)
%computes IRF pattern
lag_vec=0:lags-1;



Sigma_lagged=exp(-((lag_vec-b).^2)./c)...
    .*a;



end

