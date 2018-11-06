function cdfvals = myPLCdf(data,a,xmin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a(a<=1)=NaN;
data(data<xmin)=[];
datamin=xmin;
cdfvals=0;

cdfvals=1-(data./datamin).^(1-a);
%a=round(a);
% if a-floor(a)==0
% %%expression the Hurwitz zeta function in terms of the polygamma (psi)
% %%function
% 
% myfun=psi(a-1,data)./(((-1).^(a-1)).*factorial(a-1)');
% myminfun=psi(a-1,datamin)./(((-1).^(a-1)).*factorial(a-1)');
% 
% cdfvals=1-myfun./myminfun;
% 
% end
end

