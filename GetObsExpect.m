function [Obs Expect] = GetObsExpect(data,mycdf,params,n_bins);
%This function takes in observed data and computes its binned histogram. It
%also computes its expected (theoretical) counts according to a
%user-provided distribution function Sheyum 08/24/11


mincount=5; %minimum number of counts in any expected bin

%%first compute the frequency of unique observations linearly-spaced
%%bins,first
nbins=(nargin<=3)*20;  %number of bins to divide the data into
if nargin>3, nbins=n_bins; end

binedges=min(data)+((max(data)-min(data))/nbins)*(0:nbins);
binedges = [-Inf binedges(2:end-1) Inf];
% maxn=floor(log(max(data))/log(2)); binedges=[min(data)
% min(data)+cumsum(2.^(0:maxn)) max(data)]; binedges(binedges>max(data))=[]

Obs=histc(data,binedges);
Obs=Obs(:);
%[Obs,binedges]=poolObs(Obs,binedges,mincount); %%make sure observation
%bins contain enough counts
%%get initial expected counts from cdf of model
trimmededges=binedges(2:end-1);

% P = num2cell(params) cdfvalues=feval(mycdf,trimmededges,P{:});
% cdfvalues=[0;cdfvalues(:);1];

[xx xxx cdfvalues] = CalcCDF(trimmededges,mycdf,params,min(data),max(data));
cdfvalues=[0;cdfvalues(:);1];

Expect=sum(Obs)*diff(cdfvalues);

% Make sure expected counts are not too small
if any(Expect<mincount)
    [Expect,Obs,binedges] = poolbins(Expect,Obs,binedges,mincount);
end

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
%plot(Obs,Expect)
%%%----------------------------%%%%

    function [Obs,edges]=poolObs(Obs,edges,emin)
        i = 1;
        j = length(Obs);
        while(i<j-1 && ...
                (   Obs(i)<emin || Obs(i+1)<emin ...
                || Obs(j)<emin || Obs(j-1)<emin))
            if Obs(i)<Obs(j)
                Obs(i+1) = Obs(i+1) + Obs(i);
                i = i+1;
            else
                Obs(j-1) = Obs(j-1) + Obs(j);
                j = j-1;
            end
        end
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

