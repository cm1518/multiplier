function [ mod_params ] = params_mod( params,setup )
%takes restricted parameter draw and adds zeros for the restricted
%parameters so that original code can be used





mod_params.intercepts=params(1:setup.size_obs);
counter=setup.size_obs+1;
poly_coefficients=params(counter:counter+setup.polynomials*setup.size_obs-1);
mod_params.poly_coefficients=reshape(poly_coefficients,setup.size_obs,setup.polynomials);
counter=counter+setup.polynomials*setup.size_obs;
mod_params.beta_diag_neg=params(counter:counter+setup.size_obs-setup.index_unrestricted);
counter=counter+setup.size_obs-setup.index_unrestricted+1;

mod_params.alpha_gen_neg=zeros(setup.size_obs,max(setup.num_gaussian));
for jj=1:max(setup.num_gaussian)
    mod_params.alpha_gen_neg(setup.num_gaussian>=jj,jj)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end





%note that beta_gen is restricted to be the same for all GBs
mod_params.beta_gen_neg=zeros(setup.size_obs,1);

    mod_params.beta_gen_neg=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;



mod_params.b_gen_neg=zeros(setup.size_obs,max(setup.num_gaussian));
for jj=1:max(setup.num_gaussian)
    mod_params.b_gen_neg(setup.num_gaussian>=jj,jj)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end



mod_params.c_gen_neg=ones(setup.size_obs,max(setup.num_gaussian)); %note that here the default value is one to avoid division by zero
for jj=1:max(setup.num_gaussian)
    mod_params.c_gen_neg(setup.num_gaussian>=jj,jj)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end



mod_params.beta_diag_pos=params(counter:counter+setup.size_obs-setup.index_unrestricted);
counter=counter+setup.size_obs-setup.index_unrestricted+1;



mod_params.alpha_gen_pos=zeros(setup.size_obs,max(setup.num_gaussian));
for jj=1:max(setup.num_gaussian)
    mod_params.alpha_gen_pos(setup.num_gaussian>=jj,jj)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end

%note that beta_gen is restricted to be the same for all GBs
mod_params.beta_gen_pos=zeros(setup.size_obs,1);

    mod_params.beta_gen_pos=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;


mod_params.b_gen_pos=zeros(setup.size_obs,max(setup.num_gaussian));
for jj=1:max(setup.num_gaussian)
    mod_params.b_gen_pos(setup.num_gaussian>=jj,jj)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end

mod_params.c_gen_pos=ones(setup.size_obs,max(setup.num_gaussian)); %note that here the default value is one to avoid division by zero
for jj=1:max(setup.num_gaussian)
    mod_params.c_gen_pos(setup.num_gaussian>=jj,jj)=params(counter:(counter+sum(setup.num_gaussian>=jj)-1));
    counter=counter+(sum(setup.num_gaussian>=jj));
    
end







end

