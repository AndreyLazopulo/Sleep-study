function [cdfvals] = myPLwithHyperExpoCdf(data,params,xmin,xmax,K)
%%Cumulative distribution function of a composite Brownian walk model made
%%up of K exponential terms. Assumes the model cdf is written as
%% cdf=1-sum(p_i*exp(-lambda_i*(x-xmin)))  for i=1 to K. Notice all p_i's
%%lambda_i's are positive in this expression.

cdfvals=[];
if isempty(xmin)|| xmin<0, xmin=min(data(data>=0)); end
if isempty(xmax), xmax=max(data(data>=0)); end

data(data<xmin|data>xmax)=[];   %%exclude data that are outside range
params(params<0)=[]; %%exclude normalization constants that are negative

if numel(params)~=2*K+1
    disp('myPLwithHyperExpoCdf:Expected number of parameters not supplied')
    disp('Parameters must be zero or positive')
    return
end

a=params(1); %the first parameter is alpha
wpdf=params(2:K+1); %weight of each term in the pdf
lam=params(K+2:end); % the lambda exponents

cdfvals=0; %%initialize
for jj=1:K
    wcdf=-wpdf(jj)*(lam(jj))^(a-1); %% rescaling wpdf to get weight of cdf term
    gamnumer=gamma_incomplete(data*lam(jj),1-a)-gamma_incomplete(xmin*lam(jj),1-a);
    newterm=wcdf*gamnumer;
    cdfvals=cdfvals+newterm;
end

end