function [W] = ageweightedhs(lambda,i,N)

W=(lambda.^(i-1).*(1-lambda))./(1-lambda.^N); %compute weight for a given date and size of sample
end

