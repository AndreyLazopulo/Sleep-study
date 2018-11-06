function [cdfvals] = myBoundExpCDF(data,mu,xmin,xmax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isempty(xmin)==1, xmin=min(data); end;
if isempty(xmax)==1, xmax=max(data); end;
data(data<xmin|data>xmax)=[];   %exclude data that outside defined range


denom=exp(-xmin/mu)-exp(-xmax/mu);
numer=exp(-xmin/mu)-exp(-data./mu);

cdfvals = numer./denom;

% Force values to 0 or 1 outside range.
cdfvals(data < xmin) = 0;  cdfvals(data > xmax) = 1.0;

end

