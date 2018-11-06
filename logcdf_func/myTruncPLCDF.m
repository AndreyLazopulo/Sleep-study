function cdfvals = myTruncPLCDF(data,a,xmin,xmax)
%Cumulative distribution function for the truncated (at both ends) Pareto
%distribution. Sheyum 08/27/11
% 
data(data<xmin|data>xmax)=[];  %exclude data outside designated range

normal=xmax^(1-a)-xmin^(1-a);
cdfvals=(data.^(1-a)-xmin^(1-a))./normal;

end

