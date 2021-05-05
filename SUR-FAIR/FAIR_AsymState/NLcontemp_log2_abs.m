function [ result ] = NLcontemp_log( x, alpha_,beta_,m )
%evaluates the log of the function x*exp(alpha_*x^2+beta_)-m
%taking into account the sign of x
if m>0
%result=x.*exp(alpha_*x.^2+beta_)-m;
result=x+alpha_*abs(exp(x))+beta_-log(m);
elseif m<0
   result=x+alpha_*abs(exp(x))+beta_-log(-m);
    
else
    disp('error')
end
end

