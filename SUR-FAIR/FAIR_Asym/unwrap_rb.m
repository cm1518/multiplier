function G=unwrap_rb(y,setup)

%order variables as cte, GB1neg, GB1pos, GB2neg etc.., cov matrix, VAR coefs
% Ordering of paramters in SUR: c_G, GB1_Gneg, GB1_Gpos,  GB_2Gneg,  GB_2Gpos, c_Y, GB1_Yneg, GB1_Ypos,  GB_2Yneg,  GB_2Ypos, cov_VAR, VAR_coef

%only works for same nbGB across Y and G
nbGB=setup.num_Gaussian{1};

counter=1;
tp=[];
tp(1)=y(1);
for counter=1:2 %correspond to the pos and neg
tp=[tp; y(2+counter:2*nbGB*2:2+counter+2*nbGB*2*2+1)];
end

tp=[tp; y(2)];
for counter=1:2
tp=[tp; y(2+2*nbGB+counter:2*nbGB*2:2*nbGB+2+counter+2*nbGB*2*2+1)];
end


%add in VAR coefs and variance
G=[tp; y(end-6:end)];

end












