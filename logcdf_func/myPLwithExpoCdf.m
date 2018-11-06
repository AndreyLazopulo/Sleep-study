function [cdfvals] = myPLwithExpoCdf(data,params,xmin,xmax,K)
%%Cumulative distribution function of a model whose pdf is represented as a product
%%of a power law,x^(-Alpha), and K exponential terms:
%%% x^(-alpha)*(weight1*exp(-lambda1*x)+weight2*exp(-lambda2*x)+...

%%% The cdf=a1*Incomple gamma func1/Incomple gamma0+
%%%             a2*Incomple gamma func2/Incomple gamma0+...

%%% IMPORTANT: In 'params', param(1) is assumed to be alpha, params(2:k+1) are assumed
%%% to be the a's and params(k+2:2*k+1) are assumed to be the lambda's

cdfvals=[];

if isempty(xmin)|| xmin<0, xmin=min(data(data>=0)); end
if isempty(xmax), xmax=max(data(data>=0)); end

data(data<xmin|data>xmax)=[];   %%exclude data that are outside range

%%%%%%Start of checks and balances
if isempty(data)
    disp('myPLwithExpoCdf: All data must be greater than zero')
    disp('myPLwithExpoCdf: No data found between xmin and xmax')
    return
end

if numel(params)~=2*K+1
    disp('myPLwithExpoCdf:Expected number of parameters not supplied')
    return
end

if any(params(K+2:end)<0)
    disp('myPLwithExpoCdf: Exponents must be positive')
    return
end
%%%%End of checks and balances

%%%%Calculate the CDF
data=sort(data); cdfvals=0;
Alpha=1-params(1);
a=params(2:K+1); lam=params(K+2:2*K+1);
for q=1:K
    z=lam(q)*data; zmin=lam(q)*xmin; zmax=lam(q)*xmax;
    gammnumer=a(q)*(gamma_incomplete(z,Alpha)-gamma_incomplete(zmin,Alpha));
    gammdenom=gamma_incomplete(zmax,Alpha)-gamma_incomplete(zmin,Alpha);
    cdfvals=cdfvals+gammnumer./gammdenom;
end

end
