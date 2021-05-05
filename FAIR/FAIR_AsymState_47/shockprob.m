function [ probs ] = shockprob( add_matrices, setup )
%calculates the probability of a positive shock and plots those
%probabilities
probs=mean(add_matrices,3);

for jj=1:setup.size_obs
figure(jj)
    plot(probs(jj,:),'Linewidth',2)
    
end
end

