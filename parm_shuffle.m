function [amod,lmod] = parm_shuffle(a,l)
[Nr,Nc]=size(a);
idx=nchoosek([1:Nc],Nc-1);
amod=[]; lmod=[];
for i=1:Nr
    amod=[amod;reshape(a(i,idx),Nc,Nc-1)];
    lmod=[lmod;reshape(l(i,[idx,flipud([1:Nc]')]),Nc,Nc)];    
end
