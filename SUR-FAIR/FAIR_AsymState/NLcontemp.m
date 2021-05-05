function [ result ] = NLcontemp( x, alpha_,beta_,m )
%evaluates the function x*exp(alpha_*x^2+beta_)-m


result=x.*exp(alpha_*x.^2+beta_)-m;

end

