function G=wrap_rb(setup)

%order variables as cte, a, b, c, cov matrix, VAR coefs
y=setup.initial_parameter;

nbGB=setup.num_Gaussian{1};
counter=1;
tp_cte=y(1);
tp_a=y(counter+1:3:counter-2+2*nbGB*3);
tp_b=y(1+counter+1:3:1+counter-2+2*nbGB*3);
tp_c=y(2+counter+1:3:2+counter-2+2*nbGB*3);


nbGB=setup.num_Gaussian{2};
counter=1;
tp2_cte=y(2+3*nbGB*2);
tp2_a=y(2+3*nbGB*2+counter:3:counter-2+1+4*nbGB*3);
tp2_b=y(1+2+3*nbGB*2+counter:3:1+counter-2+1+4*nbGB*3);
tp2_c=y(2+2+3*nbGB*2+counter:3:2+counter-2+1+4*nbGB*3);

G=[tp_cte; tp2_cte; tp_a; tp2_a; tp_b; tp2_b; tp_c; tp2_c];

%add in VAR coefs and variance
G=[G;y(end-6:end)];

end












