function [ A ] = VAR_MA_implied( setup,j )
%calculates the MA coefficient of order j implied by a first order VAR
%(coefficients for which are part of setup)

A=setup.VARsymA^(j)*setup.VARsymchol;
end

