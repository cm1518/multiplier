function rigid=rigidfunc(Unemp,Unemp_n)
global gamma0 gamma1 gamma2 U_ss do_constantelast

if do_constantelast
    rigid = (1-Unemp).^gamma2 ./ (1-U_ss).^gamma2;
else
    rigid = gamma0+ 1/gamma2*(gamma1*Unemp).^(-gamma2);
end