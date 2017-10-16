function [gof pval] = Gtest(obsvals,expctdvals,noparams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gof=[]; pval=[];

size_obs=size(obsvals);
size_exp=size(expctdvals);

if size_obs~=size_exp;
    expctdvals=expctdvals';
end

gof=2*sum(obsvals.*log(obsvals./expctdvals));
degfreedom=numel(obsvals)-noparams-1;

%%William's correction;
corrfac=1+(numel(obsvals)^2-1)/(6*sum(obsvals)*degfreedom);
gof=gof/corrfac;

%large p val indicates the two are from the same distribution
if degfreedom>0
pval=1-chi2cdf(gof,degfreedom);
end
end

