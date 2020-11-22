function [failurerate,LR,result,LRCCI,CC] = failurerate(Returns,VaR,confidencelvlVaR,confidencelvltest)
%This function allows to compute the faillure rate percentage of case exceed
%the VaR and compute the kupiec, we also run a chrtisoffersen test for independancy and the
%joint test of both test. 
%   Return is the vector of size Tx1 of the return want to compute the
%   faillure rate
%   VaR is the vector of same size Tx1 as a threshold, VaR must be positive
%   function return the failure rate
%   HO is failure rate = 1-confidencelvl, we dont want to reject

%kupiec
Return=-Returns; %loss function 
var=VaR;
p=1-confidencelvlVaR;
T=size(Return,1);
Return(Return>var)=1; %if higher than the var it is a fail so it equals 1
Return(Return~=1)=0; %elsewhere put 0
N=sum(Return); %total number of fail
failurerate=N/T; %failure rate

LR=2*(log((failurerate).^N * (1-failurerate).^(T-N))-log(p.^N * (1-p).^(T-N))); %kupiec test as in the paper

if isnan(LR) %ensure no nan
    LR=inf;
end

criticalvalue= chi2inv(confidencelvltest,1);
if LR>criticalvalue %result of the test
    result=1;
else
    result=0;
end

%Christoffersen
ff=0;fs=0;sf=0;ss=0; %initialize value
for i=2:T
    if Return(i,1)==1 && Return(i-1,1)==1 %fail fail condition 
        ff=ff+1;
    elseif Return(i,1)==1 && Return(i-1,1)==0
        fs=fs+1;
    elseif Return(i,1)==0 && Return(i-1,1)==1
        sf=sf+1;
    elseif Return(i,1)==0 && Return(i-1,1)==0
        ss=ss+1;
    end
end
pi0=sf./(ss+sf); %compute the different proba
pi1=ff./(fs+ff);
pi2=(sf+ff)./(ss+ff+sf+fs);
upper=((1-pi2).^(ss+fs))*(pi2.^(sf+ff));
lower=((1-pi0).^(ss))*(pi0.^(sf))*((1-pi1).^(fs))*(pi1.^(ff));
LRCCI=-2*log(upper/lower); %test


CC=LR+LRCCI; %christofferen test

end

