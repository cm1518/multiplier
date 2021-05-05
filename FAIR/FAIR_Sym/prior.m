function [ lnp ] = prior( param,setup )
%computes log prior for normal and gamma prior distributions
%i'm assuming here that the priors on the initial response (and all other parameters!) are normal... If some parameters do not have a Gaussian prior this code will be off 
[ param ] = inv_transform( param,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
[ mod_params ] = params_mod( param,setup );

%calculate responses on impact of GMA would hold on impact
initial_neg=exp(-((mod_params.b_gen_neg).^2)./mod_params.c_gen_neg).*(mod_params.beta_gen_neg);
initial_neg=initial_neg(setup.index_unrestricted:end);


temp_means=setup.normal_prior_means;
temp_means(setup.initial_negative_index)=initial_neg;


if ~isempty(setup.index_normal)
log_normal=sum(log(normpdf(param(setup.index_normal),temp_means,setup.normal_prior_std)));
else
log_normal=0;
end


if ~isempty(setup.index_gamma)
log_gamma=sum(log(pdf('Gamma',param(setup.index_gamma),setup.gamma_prior_shape,setup.gamma_prior_scale)));
else
log_gamma=0;
end


lnp=log_normal+log_gamma;


end

