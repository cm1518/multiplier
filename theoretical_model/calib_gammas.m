function resid=calib_gammas(gammas)
global gamma0 gamma1 gamma2 U_ss pibar
gamma0 = gammas(1);
gamma1 = gammas(2);
gamma2 = gammas(3);

resid(1) = 1^(1/4) - rigidfunc(U_ss,U_ss); % wages fully downwardy rigid at steady state
resid(2) = 0.97^(1/4) - rigidfunc(U_ss+0.01,U_ss); 
resid(3) = 0.96^(1/4) - rigidfunc(U_ss+0.04,U_ss); 

% resid(1) = 1^(1/4) - rigidfunc(U_ss,U_ss); % wages fully downwardy rigid at steady state
% resid(2) = 0.97^(1/4) - rigidfunc(U_ss+0.02,U_ss); 
% resid(3) = 0.96^(1/4) - rigidfunc(U_ss+0.04,U_ss); 