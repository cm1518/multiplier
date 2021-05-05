function [ acc_rate ,param_draws,log_posteriors,x_draws] = draws_standard_RW_for_Adaptive( posterior_draw,cov_matrix,scaling,old_posterior,setup,data  )
%returns posterior draws (NOT TRANSFORMED BACK TO STANDARD PARAMETER SPACE) for the standard MH algorithm with a RW propsal
%density

acceptances=0;
param_draws=zeros(setup.length_param_vector,setup.burn_in/setup.keep_draw);
%temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );   
temp_draw=posterior_draw;
param_draws(:,1)=temp_draw;
log_posteriors=zeros(1,setup.burn_in/setup.keep_draw);

x_draws=zeros(setup.state_size,setup.sample_size+1,setup.burn_in/setup.keep_draw);

    
%    if setup.likelihood==1
%       % ll_function=SSKF_wrap( params,setup,data );
%        posterior=@(params,setup,data) (prior(params,setup)+SSKF_wrap( params,setup,data ));
%        
%    elseif setup.likelihood==2
%        %ll_function=KF_wrap;
%         posterior=@(params,setup,data) (prior(params,setup)+KF_wrap( params,setup,data ));
%      
%    else
%   
%    end

[trash xdraw trash2]=posteriorwithstates(posterior_draw,setup,data);

for nn=1:setup.burn_in
%    if (nn/setup.disp_iter)==floor(nn/setup.disp_iter)
%     nn
%    end
param_prop=proposal_draw(posterior_draw,cov_matrix, scaling,setup);
param_prop=param_prop';
[post_prop x_prop trash2]=posteriorwithstates(param_prop,setup,data);
%post_prop=posterior(param_prop,setup,data);
alpha=min(1,exp(post_prop-old_posterior));
if rand<alpha
   posterior_draw=param_prop;
   old_posterior=post_prop;
   acceptances=acceptances+1;
   %temp_draw=inv_transform( posterior_draw,setup.index_log, setup.index_logit,setup.index_logit_general, setup.length_log,setup.length_logit,setup.length_logit_general,setup.logit_general_lb,setup.logit_general_ub );
   temp_draw=posterior_draw;
xdraw=x_prop;
end

if ((nn)/setup.keep_draw)==floor(((nn)/setup.keep_draw))
param_draws(:,(nn)/setup.keep_draw)=temp_draw;
log_posteriors((nn)/setup.keep_draw)=old_posterior;
x_draws(:,:,(nn)/setup.keep_draw)=xdraw;
end

end


acc_rate=acceptances/setup.burn_in;

end
