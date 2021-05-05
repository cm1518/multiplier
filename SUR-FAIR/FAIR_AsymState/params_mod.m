function [ mod_params ] = params_mod( params,setup )
%takes restricted parameter draw and adds zeros for the restricted
%parameters so that original code can be used
length_res=length(setup.index_restricted);
length_unres=setup.size_obs-length_res;

intercepts=params(1:setup.size_obs);
counter=setup.size_obs+1;
beta_diag_neg=params(counter:counter+length_unres-1);
counter=counter+length_unres;
alpha_gen_neg=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
beta_gen_neg=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
b_gen_neg=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
c_gen_neg=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
gamma_gen_neg=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;



beta_diag_pos=params(counter:counter+length_unres-1);
counter=counter+length_unres;
alpha_gen_pos=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
beta_gen_pos=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
b_gen_pos=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
c_gen_pos=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
gamma_gen_pos=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;
c_state_neg=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;

c_state_pos=params(counter:counter+setup.size_obs*length_unres-1);
counter=counter+setup.size_obs*length_unres;

square_inds=zeros(setup.size_obs*length_unres,1);

for jj=1:length_unres
    
    square_inds(1+(jj-1)*setup.size_obs:(jj-1)*setup.size_obs+setup.size_obs)=[setup.size_obs*setup.index_unrestricted(jj)-setup.size_obs+1:setup.size_obs*setup.index_unrestricted(jj)]';
end

beta_diag_new_neg=zeros(setup.size_obs,1);
beta_diag_new_neg(setup.index_unrestricted)=beta_diag_neg;
alpha_gen_new_neg=zeros(setup.size_obs^2,1);
alpha_gen_new_neg(square_inds)=alpha_gen_neg;
beta_gen_new_neg=zeros(setup.size_obs^2,1);
beta_gen_new_neg(square_inds)=beta_gen_neg;
b_gen_new_neg=zeros(setup.size_obs^2,1);
b_gen_new_neg(square_inds)=b_gen_neg;
c_gen_new_neg=zeros(setup.size_obs^2,1);
c_gen_new_neg(square_inds)=c_gen_neg;
gamma_gen_new_neg=zeros(setup.size_obs^2,1);
gamma_gen_new_neg(square_inds)=gamma_gen_neg;


beta_diag_new_pos=zeros(setup.size_obs,1);
beta_diag_new_pos(setup.index_unrestricted)=beta_diag_pos;
alpha_gen_new_pos=zeros(setup.size_obs^2,1);
alpha_gen_new_pos(square_inds)=alpha_gen_pos;
beta_gen_new_pos=zeros(setup.size_obs^2,1);
beta_gen_new_pos(square_inds)=beta_gen_pos;
b_gen_new_pos=zeros(setup.size_obs^2,1);
b_gen_new_pos(square_inds)=b_gen_pos;
c_gen_new_pos=zeros(setup.size_obs^2,1);
c_gen_new_pos(square_inds)=c_gen_pos;
gamma_gen_new_pos=zeros(setup.size_obs^2,1);
gamma_gen_new_pos(square_inds)=gamma_gen_pos;

c_state_new_neg=zeros(setup.size_obs^2,1);
c_state_new_neg(square_inds)=c_state_neg;
c_state_new_pos=zeros(setup.size_obs^2,1);
c_state_new_pos(square_inds)=c_state_pos;

mod_params=[intercepts;beta_diag_new_neg;alpha_gen_new_neg;beta_gen_new_neg;b_gen_new_neg;c_gen_new_neg;gamma_gen_new_neg;...
    beta_diag_new_pos;alpha_gen_new_pos;beta_gen_new_pos;b_gen_new_pos;c_gen_new_pos;gamma_gen_new_pos;c_state_new_neg;c_state_new_pos];

end

