function [cdfvals] = myExpoCdf(data,mu,xmin)
%This function calculates the cumulative distribution for an exponential
%function. 
%Sheyum 08/24/11

data(data<xmin)=[];   %exclude data that is less than the assigned min val.
mu(mu <= 0) = NaN;

datamin=xmin;
x=(data-datamin)./mu;
cdfvals = 1 - exp(-x);

% Force a zero for all negative values of the argument.
cdfvals(x < 0) = 0;

end

