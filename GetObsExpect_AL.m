function [Obs Expect] = GetObsExpect_AL(data,mycdf,params,min_bin)
%This function takes in observed data and computes its binned histogram. It
%also computes its expected (theoretical) counts according to a
%user-provided distribution function Sheyum 08/24/11


mincount=5; %minimum number of counts in any expected bin

%%first compute the frequency of unique observations linearly-spaced
%%bins,first
% nbins=(nargin<=3)*20;  %number of bins to divide the data into
% if nargin>3, nbins=n_bins; end

[x y cdfvalues] = CalcCDF_forKS(data,mycdf, params,min(data),max(data));
% binedges=min(data)+((max(data)-min(data))/nbins)*(0:nbins);
% binedges = [-Inf binedges(2:end-1) Inf];
% x=binning_function(data,'constant',min_bin);
% maxn=floor(log(max(data))/log(2)); binedges=[min(data)
% min(data)+cumsum(2.^(0:maxn)) max(data)]; binedges(binedges>max(data))=[]

[Obs,edges]=histcounts(data,x);
% figure; histogram(data,x);
% figure; loglog((edges(1:end-1)+edges(2:end))./2,Obs)
Obs=Obs(:);
[Obs,x]=poolObs(Obs,x,10); %%make sure observation
% figure;histogram(data,x);
[Obs,x]=histcounts(data,x);
% loglog(binedges(1:end-1),Obs)
%bins contain enough counts
%%get initial expected counts from cdf of model
% trimmededges=binedges(2:end-1);
trimmededges=x;

% P = num2cell(params) cdfvalues=feval(mycdf,trimmededges,P{:});
% cdfvalues=[0;cdfvalues(:);1];

[xx xxx cdfvalues] = CalcCDF(trimmededges,mycdf,params,min(data),max(data));
% cdfvalues=[0;cdfvalues(:);1];
% if max(cdfvalues)<0.9999
%     cdfvalues=[cdfvalues(:);0.9999];
% end

Expect=sum(Obs)*diff(cdfvalues);

% Make sure expected counts are not too small
% if any(Expect<mincount)
%     [Expect,Obs,binedges] = poolbins(Expect,Obs,x,mincount);
% end
% figure; histogram(data,binedges);
%%if there is still discrepancy in size between Obs and Expect, then
%%shorten the larger one
if length(Obs)>length(Expect)
    differ=length(Obs)-length(Expect);
    Obs=Obs(1:end-differ);
end
if length(Obs)<length(Expect)
    differ=-length(Obs)+length(Expect);
    Expect=Expect(1:end-differ);
end
% figure;
% plot(Obs,Expect,'o',Obs,Obs)
%%%----------------------------%%%%

function [Obs,edges]=poolObs(Obs,edges,emin)
        i = 1;
        j = length(Obs);
        while i<j-1
            if Obs(i)<5
                S=0; kk=0;
                while S<emin && i+kk<j
                    kk=kk+1;
                    S=sum(Obs(i:i+kk));
                end
                edges(i+1:i+kk)=NaN;
                i = i+kk+1;
            else
                i=i+1;
            end
        end
        edges=edges(~isnan(edges));
    end
%%THIS FUNCTION IS STRAIGHT FROM MATLAB'S CHI2GOF FUNCTION
    function [Exp,Obs,edges] = poolbins(Exp,Obs,edges,emin)
        %POOLBINS Check that expected bin counts are not too small
        
        % Pool the smallest bin each time, working from the end, but avoid
        % pooling everything into one bin.  We will never pool bins except
        % at either edge (no two internal bins will get pooled together).
        i = 1;
        j = length(Exp);
        while(i<j-1 && ...
                (   Exp(i)<emin || Exp(i+1)<emin ...
                || Exp(j)<emin || Exp(j-1)<emin))
            if Exp(i)<Exp(j)
                Exp(i+1) = Exp(i+1) + Exp(i);
                Obs(i+1) = Obs(i+1) + Obs(i);
                i = i+1;
            else
                Exp(j-1) = Exp(j-1) + Exp(j);
                Obs(j-1) = Obs(j-1) + Obs(j);
                j = j-1;
            end
        end
        
        % Retain only the pooled bins
        Exp = Exp(i:j);
        Obs = Obs(i:j);
        edges(j+1:end-1) = [];  % note j is a bin number, not an edge number
        edges(2:i) = [];        % same for i
        
        % Warn if some remaining bins have expected counts too low
        if any(Exp<emin)
            warning('stats:chi2gof:LowCounts',...
                ['After pooling, some bins still have low expected counts.\n'...
                'The chi-square approximation may not be accurate']);
        end
        
    end
%%%----------------------------%%%%
end

