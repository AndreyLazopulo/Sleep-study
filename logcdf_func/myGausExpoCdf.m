function [cdfvals] = myGausExpoCdf(data,params,xmin,xmax,K)

cdfvals=[];
if isempty(xmin)|| xmin<0, xmin=min(data(data>=0)); end
if isempty(xmax), xmax=max(data(data>=0)); end

data(data<xmin|data>xmax)=[];   %%exclude data that are outside range
% params(params(1:K)<0)=[]; %%exclude normalization constants that are negative 


if numel(params)~=2*K+1
    disp('myGausExpoCdf:Expected number of parameters not supplied') 
    disp('Parameters must be zero or positive')
    return
end

z=data; 
zmax=xmax;
zmin=xmin;
p=params(1:K);
mu=params(K+1);
sig=params(K+2);
lam=params(K+3:end);

cdfvals=p(1)*(erf((mu-z)/(sqrt(2)*sig))-erf((mu-zmin)/(sqrt(2)*sig)));  %%initialize with first term
cdfvals=cdfvals./(erf((mu-zmax)/(sqrt(2)*sig))-erf((mu-zmin)/(sqrt(2)*sig)));
for ii=1:K-1
    const=exp(-lam(ii)*zmax)-exp(-lam(ii)*zmin);
    newterm=p(ii+1)*(exp(-lam(ii)*z)-exp(-lam(ii)*zmin))./const;
    cdfvals=cdfvals+newterm;
end
    
end

