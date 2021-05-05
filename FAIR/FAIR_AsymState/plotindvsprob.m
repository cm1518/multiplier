load eps_file
ind=epsilon>0;
[ probs ] = shockprob( add_matrices, setup );
tend=50;
figure;
subplot(2,1,1)
plot(1:tend,ind(1,1:tend))
subplot(2,1,2)
plot(1:tend,probs(1,1:tend))