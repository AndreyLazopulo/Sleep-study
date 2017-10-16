function [xvals yvals cdfvals] = CalcCDF_forKS(empdistr,model, params,xmin,xmax);
%%%The function calculates the cumulative distribution function (cdf) of a
%%%given data set and the theoretical cdf using given params and name of
%%%model. The computed empirical and theoretical cdfs are returned.
%%%%Sheyum 12/20/11.

%%%For the KS test. CDF is computed only on unique values of x Sheyum
%%%3/3/2012

%%%some matlab function
[yCDF,xCDF,n,emsg,eid] = cdfcalc(empdistr);
xvals=xCDF; yvals=yCDF(1:end-1);
% xvals=[xCDF; xCDF(end)+1]; yvals=yCDF;
unqx=xvals;

%[yvals,xvals]=ecdf(empdistr);


%compute the cumulative distr of the data
% x=reshape(sort(empdistr),numel(empdistr),1); unqx=unique(x); temp =
% hist(x,unqx)'./length(x); distr = [[unqx; unqx(end)+1] 1-[0;
% cumsum(temp)]]; distr(distr(:,2)<10^-10,:) = []; unqx=distr(:,1);
% xvals=unqx;yvals=1-distr(:,2); %%%%*************
% x=reshape(sort(empdistr),numel(empdistr),1); x=x(x>=xmin); bins=[-Inf; x;
% Inf]; counts=histc(x,bins); temp=cumsum(counts)./sum(counts); distr =
% [[x; x(end)+1] 1-temp(1:end-1)]; distr(distr(:,2)<10^-10,:) = [];
% xvals=distr(:,1);yvals=1-distr(:,2);

%%calculate the correct cumulative distr function
switch model,
    case 'Exponential',
        cdfvals = myExpoCdf(xvals,params(1),xmin);
    case 'Expo_bound',
        cdfvals = myBoundExpCDF(xvals,params(1),xmin,xmax);
    case '2Exponential',
        cdfvals = myHyperExpoCdf(xvals,params,xmin,xmax,2);
    case '3Exponential',
        cdfvals = myHyperExpoCdf(xvals,params,xmin,xmax,3);
    case '4Exponential',
        cdfvals = myHyperExpoCdf(xvals,params,xmin,xmax,4);
    case '5Exponential',
        cdfvals = myHyperExpoCdf(xvals,params,xmin,xmax,5);
    case '6Exponential',
        cdfvals = myHyperExpoCdf(xvals,params,xmin,xmax,6);
    case 'PL2Exponential'
        cdfvals = myPLwithExpoCdf(xvals,params,xmin,xmax,2);
    case 'PL3Exponential'
        cdfvals = myPLwithExpoCdf(xvals,params,xmin,xmax,3);
    case 'PL4Exponential'
        cdfvals = myPLwithExpoCdf(xvals,params,xmin,xmax,4);
    case 'PL5Exponential'
        cdfvals = myPLwithExpoCdf(xvals,params,xmin,xmax,5);
    case 'PL6Exponential'
        cdfvals = myPLwithExpoCdf(xvals,params,xmin,xmax,6);
    case 'StrExpo'
        cdfvals = myStrExpoCdf(xvals,params,xmin,xmax,1);
    case 'Str1Exponential'
        cdfvals = myStrExpoCdf(xvals,params,xmin,xmax,2);
    case 'Str2Exponential'
        cdfvals = myStrExpoCdf(xvals,params,xmin,xmax,3);
    case 'Str3Exponential'
        cdfvals = myStrExpoCdf(xvals,params,xmin,xmax,4);
    case 'Str4Exponential'
        cdfvals = myStrExpoCdf(xvals,params,xmin,xmax,5);
    case 'Str5Exponential'
        cdfvals = myStrExpoCdf(xvals,params,xmin,xmax,6);
    case 'PLplus2Exponential'
        cdfvals = myPLplusHyperExpoCdf(xvals,params,xmin,xmax,2);
    case 'PLplus3Exponential'
        cdfvals = myPLplusHyperExpoCdf(xvals,params,xmin,xmax,3);
    case 'PLplus4Exponential'
        cdfvals = myPLplusHyperExpoCdf(xvals,params,xmin,xmax,4);
    case 'PLplus5Exponential'
        cdfvals = myPLplusHyperExpoCdf(xvals,params,xmin,xmax,5);
    case 'PLplus6Exponential'
        cdfvals = myPLplusHyperExpoCdf(xvals,params,xmin,xmax,6);
    case 'Log-Norm',
        cdfvals = myLogNormCDF(xvals,params(1),params(2),xmin,xmax);
    case 'Gamma',
        cdfvals = myGammaCDF(xvals,params(1),params(2),xmin,xmax);
    case 'PowerLaw_bound',
        cdfvals = myTruncPLCDF(xvals,params(1),xmin,xmax);
    case 'PLExpoCutff',
        cdfvals = myPLExpoCutCDF(xvals,params(1),params(2),xmin,xmax);
    case 'genLogNorm'
        cdfvals = mygenLogNormCDF(xvals,params(1),params(2),params(3),xmin,xmax);
    case 'genWeibull'
        cdfvals = mygenWeibullCDF(xvals,params(1),params(2),params(3),xmin,xmax);
    case 'PLLog'
        cdfvals=myPLlogCutCDF(xvals,params(1),params(2),params(3),xmin,xmax);
    case 'PLExpoLog'
        cdfvals=myPLExpologCDF(xvals,params(1),params(2),params(3),params(4),xmin,xmax);
    case 'PowerLaw',
        cdfvals = myPLCdf(xvals,params(1),xmin);
    otherwise,
        disp('Model name not recognized. Returning empty function')
        cdfvals=[];
        return
end

end

