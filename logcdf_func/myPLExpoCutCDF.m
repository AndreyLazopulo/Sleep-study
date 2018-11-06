function cdfvalus = myPLExpoCutCDF(distr,a,lam,xmin,xmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isempty(xmin)==1, xmin=min(distr), end;
if isempty(xmax)==1, xmax=max(distr), end;
distr(distr<xmin|distr>xmax)=[];

x=distr;

up=gamma_incomplete(xmin*lam,1-a)-gamma_incomplete(x*lam,1-a);
down=gamma_incomplete(xmin*lam,1-a)-gamma_incomplete(xmax*lam,1-a);
cdfvalus=up./down;
cdfvalus=reshape(cdfvalus,length(cdfvalus),1);

end

