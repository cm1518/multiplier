
close all
clear all



load sym_paper

setup.number_of_draws=setup.keep_draw*size(draws,2);
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,101,1);

indicator=0;
index=0;
IRFs=[];
hor=20;
upper=95;
lower=5;

for size_of_shock=[-1 1]
    index=index+1;
    for kk=1:length(indices_for_draws)
        
      
        shock_contemp=zeros(setup.size_obs,1);
        shock_contemp(setup.index_unrestricted)=size_of_shock;
        ind_vec=[zeros(1,setup.lags)];
        for jj=1:setup.lags+1
            epsilon=[zeros(setup.size_obs,jj-1) shock_contemp zeros(setup.size_obs,setup.lags+1-jj)];
            [ Sigma, intercept] = unwrap_NL_IRF( draws(:,indices_for_draws(kk)),epsilon,setup ,ind_vec);
            if jj==1
                ind_vec=[ind_vec(1,1:end-1) indicator];
            else
                ind_vec=[ind_vec(1,2:end) 0];
            end
            
            IRFs(:,jj,kk)=Sigma(:,setup.index_unrestricted,jj);
            IRFsgn(:,jj,kk,index)=Sigma(:,setup.index_unrestricted,jj);
            
        end
        mult(:,kk)=sum(IRFs(4,1:hor,kk))./sum(IRFs(2,1:hor,kk));
        clear ind_vec
        
        
        
    end
    
    IRFsmed(:,:,index)=squeeze(median(IRFs,3));
    IRFslow(:,:,index)=squeeze(prctile(IRFs,lower,3));
    IRFshigh(:,:,index)=squeeze(prctile(IRFs,upper,3));
    mult_med(:,:)=squeeze(median(mult,2));
    mult_lower(:,:)=squeeze(prctile(mult,lower,2));
    mult_upper(:,:)=squeeze(prctile(mult,upper,2));
    
    
end




display('multiplier: 5%, median, 95%')
[mult_lower mult_med mult_upper]

results_sym.lower=mult_lower;
results_sym.med=mult_med;
results_sym.upper=mult_upper;

save sym_results_for_table results_sym
