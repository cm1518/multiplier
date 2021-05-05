function [ A B C D ] = wrapper_file( param,setup,data )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
A=eye(2);
B=zeros(2,2);
C=[param(1) 0; 0 param(2)];
D=eye(2);
end

