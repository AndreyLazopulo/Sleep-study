function cdfvalus = myTruncPLExpoCdf(data,a,lam,xmin,xmax)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if isempty(xmin)==1, xmin=min(data); end
if isempty(xmax)==1, xmax=max(data); end
data(data<xmin|data>xmax)=[];
x=data;

up=gammainc(xmin*lam,1-a,'upper')-gammainc(x*lam,1-a,'upper');
down=gammainc(xmin*lam,1-a,'upper')-gammainc(xmax*lam,1-a,'upper');
cdfvalus=up./down;




end

