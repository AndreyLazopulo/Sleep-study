function [D KSp] = KolmogSmirnov(empdistr,theordistr, params,xmin,xmax)


[x y cdfvals] = CalcCDF_forKS(empdistr,theordistr, params,xmin,xmax);
n=length(y); j=(1:101)';   %%%%these are used inside the looop


%%plot the empirical and theoretical ccdfs
x=reshape(x,n,1); y=reshape(y,n,1); cdfvals=reshape(cdfvals,length(cdfvals),1);
distr(:,1)=x; distr(:,2)=1-y; ccdfvals = [distr(:,1) 1-cdfvals];

%handl1=LogPlotDistr(distr,ccdfvals);

%%%initialize a few variables
D=[]; KSp=[];

if ~isempty(cdfvals)
    %%%%calculate the KS statistic for 2-sided test
    y=reshape(y,n,1); cdfvals=reshape(cdfvals,length(cdfvals),1);
    D=max(abs(cdfvals - y));
    
    %%% Use the asymptotic Q-function to approximate the 2-sided
    %%% P-value, the probablity that the empirical cdf is from the
    %%% theoretical cdf. The computation is according to Matlab's
    %%% kstest2.m, which is essentially outlined in Numerical Recipes
    %%% Chaps 6.14.12 and 14.3.3
    lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * D , 0);
    KSp  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
    KSp =  min(max(KSp, 0), 1);
else
    D=[]; KSp=[];
    
end



end

