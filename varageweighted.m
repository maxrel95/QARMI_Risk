function [var] = varageweighted(data,lambda,varlvl)
%this function compute the age weighted var
%   DATA        is a matrix Tx4 of 4 columns first one need to be the age of the data,
%               the second one is the return, and 2 columns of zeros, in 3rd one
%               we compute the weight and in the 4th we compute the cumulative sum
%  lambda       muste bet between]0,1[ creator based at 0.98
%  varlvl       is a row vector, with the var lvl we consider 5%,1%,...

N=size(varlvl,2);
T=length(data);
var=zeros(1,N);
for i=1:N
    test2=sortrows(data,2); %sort the data according to the second column be keeps the corresponding value other col
    test2(:,3)=ageweightedhs(lambda,test2(:,1),T); %get the weight 
    test2(:,4)=cumsum(test2(:,3)); 
    var(1,i)=-test2(find( test2(:,4) > varlvl(1,i), 1 ),2); %find when we exceed the varlvl in the cum sum and return the
end %corresponding VaR.
end

