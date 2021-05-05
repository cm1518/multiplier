% load test_SUR
% % draws=setup.initial_parameter;
% draws=draws(:,1:7001);

% draws=draws(:,501:end);

setup.number_of_draws=setup.keep_draw*size(draws,2);
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,100,1);

YGr=1; %ATTENTION, here it's 1, because things have been ex-ante normalized

%indicator:
ind=-2;

hor=20;

IRF_neg=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));
IRF_pos=zeros(setup.dep_variables,setup.lags_shocks,length(indices_for_draws));

Mmax=[]; Msum=[];

for kk=1:length(indices_for_draws)
    params=draws(:,indices_for_draws(kk));
    counter=1;
    for jj=1:setup.dep_variables
        
        %unwrapping parameters
        
        intercept{jj}=params(counter);
        counter=1+counter;
        
        
        
        lagged_y_param{jj}=params(counter:counter+setup.lags_endo-1);
        counter=counter+setup.lags_endo;
        
        lagged_exo_param{jj}=params(counter:counter+setup.lags_exo*setup.dim_exo-1); %I assume that parameters here are ordered by variable first, then by lag
        counter=counter+setup.lags_exo*setup.dim_exo;
        
        
        if jj==1
            counter_linear=counter; %number of 'linear' (except AR coefficients) plus 1, same for all dep. variables
        end
        
        NL_a_neg{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_b_neg{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_c_neg{jj}=zeros(setup.num_Gaussian{jj},1);
        
        
        
        NL_a_pos{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_b_pos{jj}=zeros(setup.num_Gaussian{jj},1);
        
        NL_c_pos{jj}=zeros(setup.num_Gaussian{jj},1);
        
        for ll=1:setup.num_Gaussian{jj}
            NL_a_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_b_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_c_neg{jj}(ll)=params(counter);
            counter=counter+1;
            NL_a_pos{jj}(ll)=params(counter);
            counter=counter+1;
            NL_b_pos{jj}(ll)=params(counter);
            counter=counter+1;
            NL_c_pos{jj}(ll)=params(counter);
            counter=counter+1;
        end
        
        
        gamma_neg{jj}=params(counter);
        counter=counter+1;
         
        beta_neg{jj}=params(counter);
        counter=counter+1;
         
        gamma_pos{jj}=params(counter);
        counter=counter+1;
        beta_pos{jj}=params(counter);
        counter=counter+1;
        
    end
    
    
    for jj=1:setup.dep_variables
        Sigma_lagged_neg{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
        for gg=1:setup.num_Gaussian{jj}
            [ Sigma_lagged_neg{jj}(gg,:) ] = NL_IRF(NL_a_neg{jj}(gg),NL_b_neg{jj}(gg),NL_c_neg{jj}(gg),gamma_neg{jj},beta_neg{jj},setup.lags_shocks,ind);
            
        end
        
        Sigma_lagged_pos{jj}=zeros(setup.num_Gaussian{jj},setup.lags_shocks);
        for gg=1:setup.num_Gaussian{jj}
            [ Sigma_lagged_pos{jj}(gg,:) ] = NL_IRF(NL_a_pos{jj}(gg),NL_b_pos{jj}(gg),NL_c_pos{jj}(gg),gamma_pos{jj},beta_pos{jj},setup.lags_shocks,ind);
            
        end
        
        
        IRF_neg(jj,:,kk)=sum(Sigma_lagged_neg{jj},1);
        IRF_pos(jj,:,kk)=sum(Sigma_lagged_pos{jj},1);
    end
    
    % Store multipliers matrix, first column with neg shocks, second column for pos shocks
    [mx_pos i_mxpos]=max(YGr*abs(IRF_pos(2,:,kk)));
    [mx_neg i_mxneg]=max(YGr*abs(IRF_neg(2,:,kk)));
    
    Mmax(kk,:)=[sign(IRF_neg(2,i_mxneg,kk)).*mx_neg./max(IRF_neg(1,:,kk)) sign(IRF_pos(2,i_mxpos,kk)).*mx_pos./max(IRF_pos(1,:,kk))];
    Msum(kk,:)=[sum(YGr*IRF_neg(2,:,kk))./sum(IRF_neg(1,:,kk)) sum(YGr*IRF_pos(2,:,kk))./sum(IRF_pos(1,:,kk))];

    
end
%moments of irfs
upper=95;
lower=5;
median_neg=prctile(IRF_neg,50,3); %first dimension is the observable, second the horizon
lower_neg=prctile(IRF_neg,lower,3);
upper_neg=prctile(IRF_neg,upper,3);

median_pos=prctile(IRF_pos,50,3);
lower_pos=prctile(IRF_pos,lower,3);
upper_pos=prctile(IRF_pos,upper,3);

mmax=squeeze(median(Mmax,1));
mmaxlow=squeeze(prctile(Mmax,lower,1));
mmaxhigh=squeeze(prctile(Mmax,upper,1));

msum=squeeze(median(Msum,1));
msumlow=squeeze(prctile(Msum,lower,1));
msumhigh=squeeze(prctile(Msum,upper,1));

figure(1),
for jj=1:setup.dep_variables
    subplot(2,2,jj)
    plot(1:setup.lags_shocks,median_pos(jj,:),'Linewidth',2,'Color','b')
    hold on, plot(1:setup.lags_shocks,lower_pos(jj,:),1:setup.lags_shocks,upper_pos(jj,:),'Color','b')
    if jj==1
        ylabel('Expansionary shock')
        title('Government spending')
    else
        title('Output')
    end
    subplot(2,2,2+jj)
    plot(1:setup.lags_shocks,median_neg(jj,:),'Linewidth',2,'Color','r')
    hold on, plot(1:setup.lags_shocks,lower_neg(jj,:),1:setup.lags_shocks,upper_neg(jj,:),'Color','r')
    if jj==1
        ylabel('Contractionary shock')
        title('(-) Government spending')
    else
        title('(-) Output')
    end
end


%Normalized IRFs
x=1:setup.lags_shocks;
figure,
for jj=1:setup.dep_variables
    
    subplot(2,2,jj), hold on
%     g2 = plot(x,0*median_pos(2,:) ,'LineWidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    error1_pos = (median_pos(jj,:)./max(median_pos(1,:)) - lower_pos(jj,:)./max(median_pos(1,:)))';
    error2_pos = (upper_pos(jj,:)./max(median_pos(1,:)) - median_pos(jj,:)./max(median_pos(1,:)))';
    H= shadedErrorBar(x,median_pos(jj,:)./max(median_pos(1,:)),[error2_pos'; error1_pos'], 'b');
    plot(median_pos(jj,:)./max(median_pos(1,:)),'Linewidth',2,'Color','b')
%     g1=plot(irfsym(jj,:)./max(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--')
    if jj==1
        ylabel('Expansionary shock')
        title('Government spending')
        legend([g1],{'Linear'}), legend('boxoff')
    else
        title('Output')
        ylim([0 1.75])
    end
    xlim([0.7 hor])
    subplot(2,2,2+jj), hold on
%     g2 = plot(x,0*median_neg(2,:) ,'LineWidth',1,'Color',[.7 .7 .7],'LineStyle','--');
    error1_neg = (median_neg(jj,:)./max(median_neg(1,:)) - lower_neg(jj,:)./max(median_neg(1,:)))';
    error2_neg = (upper_neg(jj,:)./max(median_neg(1,:)) - median_neg(jj,:)./max(median_neg(1,:)))';
    H= shadedErrorBar(x,median_neg(jj,:)./max(median_neg(1,:)),[error2_neg'; error1_neg'], 'r');
    plot(median_neg(jj,:)./max(median_neg(1,:)),'Linewidth',2,'Color','r')
%     g1=plot(irfsym(jj,:)./max(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--');
    
    %     plot(1:setup.lags_shocks,median_neg(jj,:)./max(median_neg(1,:)),'Linewidth',2,'Color','r')
    %     hold on, plot(1:setup.lags_shocks,lower_neg(jj,:)./max(median_neg(1,:)),1:setup.lags_shocks,upper_neg(jj,:)./max(median_neg(1,:)),'Color','r')
    %     plot(irfsym(jj,:)./max(irfsym(1,:)),'Linewidth',2,'Color','k','LineStyle','--')
    if jj==1
        ylabel('Contractionary shock')
        title('(-) Government spending')
    else
        title('(-) Output')
        ylim([-.0 1.75])
    end
    xlim([0.7 hor])
end


%"max" multiplier, comparing with linear GMA estimate
display('  ')
display('  ')
display('                   Neg          Pos')
display(['multiplier mmax  ', num2str([ mmax])])
display(['5% lower band:   ', num2str([ mmaxlow])]);
display(['95% upper band:  ', num2str([ mmaxhigh])]);
display(['Std dev:         ', num2str([ std(Mmax) ])]);

%"sum" multiplier, comparing with linear GMA estimate
display('  ')
display('  ')
display('                  Neg          Pos')
display(['multiplier msum  ', num2str([ msum])])
display(['5% lower band:   ', num2str([ msumlow])]);
display(['95% upper band:  ', num2str([ msumhigh])]);
display(['Std dev:         ', num2str([ std(Msum) ])]);

%
