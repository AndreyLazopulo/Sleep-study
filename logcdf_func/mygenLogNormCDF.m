function cdfvals = mygenLogNormCDF(data,alph,mu,sig,xmin,xmax )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

cdfvals=[];
data(data<xmin|data>xmax)=[];   %exclude data that are outside bounds

%%first calculate the normalization constant
zAbs=@(z)(abs((log(z)-mu)./sig)).^alph ;
f=@(z)z.^(-1).*exp(-zAbs(z)./alph);
Const=quadgk(f,xmin,xmax);
Const=1./Const;  %%%norm const determined

[unqx,numunqx] = count_unique(data);
cdfvals=zeros(length(data),1); %%%initialize
k=1;  %%initialize counter
for i=1:length(unqx)
    cdfvals(k:k+numunqx(i)-1,1)=quadgk(f,xmin,unqx(i));
    k=k+numunqx(i);
end

cdfvals=cdfvals*Const;  %%multiply cdf by normalization const

end

