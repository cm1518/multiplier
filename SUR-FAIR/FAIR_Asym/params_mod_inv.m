function [ params] = params_mod_inv( params,setup )
%takes a parameter vector from the unrestricted case and returns the
%correpsonding vector from the restricted case

%unwrapping parameter vector
intercepts=params(1:setup.size_obs);
counter=setup.size_obs+1;


beta_diag_neg=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

alpha_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

b_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

c_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

gamma_gen_neg=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_diag_pos=params(counter:counter+setup.size_obs-1);
counter=counter+setup.size_obs;

alpha_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

beta_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

b_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

c_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;

gamma_gen_pos=params(counter:counter+setup.size_obs^2-1);
counter=counter+setup.size_obs^2;





beta_diag_neg=beta_diag_neg(setup.index_unrestricted);

alpha_gen_neg=reshape(alpha_gen_neg,setup.size_obs,setup.size_obs);
alpha_gen_neg=alpha_gen_neg(:,setup.index_unrestricted);
alpha_gen_neg=alpha_gen_neg(:);

beta_gen_neg=reshape(beta_gen_neg,setup.size_obs,setup.size_obs);
beta_gen_neg=beta_gen_neg(:,setup.index_unrestricted);
beta_gen_neg=beta_gen_neg(:);

b_gen_neg=reshape(b_gen_neg,setup.size_obs,setup.size_obs);
b_gen_neg=b_gen_neg(:,setup.index_unrestricted);
b_gen_neg=b_gen_neg(:);

c_gen_neg=reshape(c_gen_neg,setup.size_obs,setup.size_obs);
c_gen_neg=c_gen_neg(:,setup.index_unrestricted);
c_gen_neg=c_gen_neg(:);

gamma_gen_neg=reshape(gamma_gen_neg,setup.size_obs,setup.size_obs);
gamma_gen_neg=gamma_gen_neg(:,setup.index_unrestricted);
gamma_gen_neg=gamma_gen_neg(:);



beta_diag_pos=beta_diag_pos(setup.index_unrestricted);

alpha_gen_pos=reshape(alpha_gen_pos,setup.size_obs,setup.size_obs);
alpha_gen_pos=alpha_gen_pos(:,setup.index_unrestricted);
alpha_gen_pos=alpha_gen_pos(:);

beta_gen_pos=reshape(beta_gen_pos,setup.size_obs,setup.size_obs);
beta_gen_pos=beta_gen_pos(:,setup.index_unrestricted);
beta_gen_pos=beta_gen_pos(:);

b_gen_pos=reshape(b_gen_pos,setup.size_obs,setup.size_obs);
b_gen_pos=b_gen_pos(:,setup.index_unrestricted);
b_gen_pos=b_gen_pos(:);

c_gen_pos=reshape(c_gen_pos,setup.size_obs,setup.size_obs);
c_gen_pos=c_gen_pos(:,setup.index_unrestricted);
c_gen_pos=c_gen_pos(:);

gamma_gen_pos=reshape(gamma_gen_pos,setup.size_obs,setup.size_obs);
gamma_gen_pos=gamma_gen_pos(:,setup.index_unrestricted);
gamma_gen_pos=gamma_gen_pos(:);





params=[intercepts;beta_diag_neg;alpha_gen_neg;beta_gen_neg;b_gen_neg;c_gen_neg;gamma_gen_neg;beta_diag_pos;alpha_gen_pos;beta_gen_pos;b_gen_pos;c_gen_pos;gamma_gen_pos];

