function [cdfvals ] = mygenWeibullCDF(data,alph,bet,lam,xmin,xmax )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cdfvals=[];
data(data<xmin|data>xmax)=[];   %exclude data that are outside bounds

nn=(alph+lam-1)/lam; z=(data./bet).^lam;
zmax=(xmax/bet)^lam; zmin=(xmin/bet)^lam; 

deno=xmin*xmax^alph*mfun('Ei',nn,zmin)-xmax*xmin^alph*mfun('Ei',nn,zmax);
numer=xmin*data.^alph*mfun('Ei',nn,zmin)-data.*mfun('Ei',nn,z)*xmin^alph;

cdfvals=(xmax^alph*data.^(-alph).*numer)/deno;
end

