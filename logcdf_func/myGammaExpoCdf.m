function [cdfvals] = myGammaExpoCdf(data,params,xmin,xmax,K)

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
alpha=1-params(K+1);
lamGamma=params(K+2);
lam=params(K+3:end);

cdfvals=p(1)*(gamma_incomplete(lamGamma*z,alpha)-gamma_incomplete(zmin*lamGamma,alpha));  %%initialize with first term
cdfvals=cdfvals./(gamma_incomplete(zmax*lamGamma,alpha)-gamma_incomplete(zmin*lamGamma,alpha));
cdfvals=cdfvals';
for ii=1:K-1
    const=exp(-lam(ii)*zmax)-exp(-lam(ii)*zmin);
    newterm=p(ii+1)*(exp(-lam(ii)*z)-exp(-lam(ii)*zmin))./const;
    cdfvals=cdfvals+newterm;
end
    
end

