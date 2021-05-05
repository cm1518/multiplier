function [ log_l indicator epsilon uvec] = likelihood_test2(data, intercept, current_matrices,lagged_matrices,eps_initial, lag_length, u )
%computes the log likelihood for the MA model
indicator=zeros(size(data,1),lag_length+size(data,2));
epsilon=zeros(size(data,1),lag_length+size(data,2));
%data is supposed to be a matrix with each column containing the
%observations at one point in time
uvec=zeros(size(data,1),size(data,2));
temp=(eps_initial>0);
indicator(:,1:lag_length)=temp;
epsilon(:,1:lag_length)=eps_initial;
%indicator equal to 0 if epsilon is negative or zero, 1 otherwise
log_l=0;
for tt=1:size(data,2)
[ log_l_temp indicator_temp epsilon_temp u_new]=likelihood_increment_u(data(:,tt), intercept,current_matrices,lagged_matrices,epsilon(:,tt:tt+lag_length-1), lag_length,indicator(:,tt:lag_length+tt-1),u(:,tt) );
log_l=log_l+log_l_temp;
indicator(:,lag_length+tt)=indicator_temp;
epsilon(:,lag_length+tt)=epsilon_temp;
uvec(:,tt)=u_new;
end
