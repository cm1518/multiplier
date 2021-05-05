function [ VAR_match ] = IRF_match( params,setup,store_responses)
%objective function for VAR responses matching below diagonal

[ Sigma_lagged ] = NL_IRF(params(1),params(2),params(3),setup.lags_shocks);
sr=store_responses;
% size(MA_h)
% size(store_responses)
% size(sr)
% size(Sigma_lagged)
VAR_match=(Sigma_lagged(:)-sr(:))'*(Sigma_lagged(:)-sr(:));
end

