function [cdfvals] = myHyperExpoCdf(data,params,xmin,xmax,K)
%%Cumulative distribution function of a composite Brownian walk model made
%%up of K exponential terms. Assumes the model cdf is written as
%% cdf=1-sum(p_i*exp(-lambda_i*(x-xmin)))  for i=1 to K. Notice all p_i's 
%%lambda_i's are positive in this expression.

cdfvals=[];
if isempty(xmin)|| xmin<0, xmin=min(data(data>=0)); end
if isempty(xmax), xmax=max(data(data>=0)); end

data(data<xmin|data>xmax)=[];   %%exclude data that are outside range
% params(params(1:K)<0)=[]; %%exclude normalization constants that are negative 


if numel(params)~=2*K
    disp('myHyperExpoCdf:Expected number of parameters not supplied') 
    disp('Parameters must be zero or positive')
    return
end

z=data-xmin; 
zmax=xmax-xmin;
p=params(1:K);
lam=params(K+1:end);

cdfvals=p(1)*(exp(-lam(1)*z)-1);  %%initialize with first term
cdfvals=cdfvals./(exp(-lam(1)*zmax)-1);
for ii=1:K-1
    const=exp(-lam(ii+1)*zmax)-1;
    newterm=p(ii+1)*(exp(-lam(ii+1)*z)-1)./const;
    cdfvals=cdfvals+newterm;
end
    
end

