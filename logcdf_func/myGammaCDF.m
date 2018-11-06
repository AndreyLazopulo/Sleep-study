function cdfvals = myGammaCDF(data,a,b,xmin,xmax)
%%The function computes the cumulative distribution function for the gamma
%%distribution. Sheyum 08/24/11
cdfvals=[];
if isempty(xmin), xmin=min(data); end
if isempty(xmax), xmax=max(data); end
data(data<xmin|data>xmax)=[];   %exclude data that is less than the assigned min val.

if xmin<0 || xmax<0
    disp('myGammaCDF:Min and Max values of the data must be positive')
    return
end

denom=gamma_incomplete(xmin/b,a)-gamma_incomplete(xmax/b,a);
cdfvals=(gamma_incomplete(xmin/b,a)-gamma_incomplete(data./b,a))./denom;
end

