function [cdfvals ] = myLogNormCDF(data,mu,sigma,xmin,xmax )
%%Function computes the cdf of a log-normal distribution defined between
%%x=xmin and x=xmax with mu and sigma

cdfvals=[];
data(data<xmin|data>xmax)=[];   %exclude data that are outside bounds
% mu(mu <= 0) = NaN; sigma(sigma <= 0) = NaN;
if sigma <= 0
    disp('myLogNormCDF:Sigma must be greater than zero');
    return
end

z=(mu-log(data))./(sqrt(2)*sigma);
zmin=(mu-log(xmin))/(sqrt(2)*sigma);
zmax=(mu-log(xmax))/(sqrt(2)*sigma);

cdfvals = (erf(zmin)-erf(z))./(erf(zmin)-erf(zmax));   %the cdf

end

