function [y,p] = qnwlognalt(ny,mu,var)
pd = makedist('Lognormal','mu',mu,'sigma',sqrt(var));
c = nodeunif(ny+1,0,1);
c = (c(1:ny)+c(2:ny+1))/2;
y = icdf(pd,c);
p = ones(ny,1)/ny;