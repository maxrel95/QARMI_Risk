function [loglik] = LLskewt(data,para)

mu=para(1); %define parameter
sigma2=para(2);
nu=para(3);
lambda=para(4);
LL=skewtloglik(data,mu,sigma2,nu,lambda); %get the loglik for each data in th sample with the corresponding parameter
loglik=-LL; %take minus because minimizer
end

