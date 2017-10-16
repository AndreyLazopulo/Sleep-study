function [R_val p_val] = LLratiotest(LL1,LL2)
%%This is the log-likelihood ratio test as described in Clauset etal SIAM
%%Review 2009. Sheyum 8/10/11

if LL2==0 | LL1==0
    p_val=1;
    R_val=0;
    return;
end

R_val=sum(LL1-LL2);   %eq. C.3 in Clauset et al

sigmaSq=mean((LL1-LL2-mean(LL1)+mean(LL2)).^2);

%%%probability of the R value being a chance result. p>0.1 indicates
%%%the difference in the LLs is not real

if R_val==0
    p_val=1; %erfc(0)
else
    p_val=erfc(abs(R_val)/sqrt(2*numel(LL1)*sigmaSq));  %eq. C.6 in Clauset
end

end

