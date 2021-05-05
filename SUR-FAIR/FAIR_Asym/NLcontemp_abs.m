function [ result ] = NLcontemp_abs( x, alpha_,beta_,m )
%evaluates the function x*exp(alpha_*x^2+beta_)-m


result=x.*exp(alpha_*abs(x)+beta_)-m;

end

