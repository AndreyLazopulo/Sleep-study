function [estAlpha AIC gof pval LL] = Binned_TruncPL(distributdata,plotting,binway)

if nargin<3
    binway=1;
elseif nargin<2
    plotting=0;
    binway=1;
end

estAlpha=[]; AIC=[]; gof=[]; pval=[];
xi=min(distributdata); xf=max(distributdata);

U=unique(distributdata);
interval=U(2)-U(1);

switch binway
    case 1
        bin_bound=[unique(distributdata);Inf];
        l2=floor(min(distributdata)):interval:max(distributdata); u2=l2+interval;
    case 2
        bin_bound=floor(min(distributdata)):interval:max(distributdata)+interval; %[unique(distributdata);Inf];
        l2=floor(min(distributdata)):interval:max(distributdata); u2=l2+interval;
end
[h,edges]=histcounts(distributdata,bin_bound);
h = reshape(h, 1, numel(h));
edges = reshape(edges, 1, numel(edges));
l=edges(1:end-1);
u=edges(2:end);
bmin=1;
specoption=optimset('Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',10000,'TolX',10^-6,'TolFun',10^-8,'Display','off','FinDiffType','central');
% indata(distributdata<xi|distributdata>xf)=[];   %exclude data outside range
estAlpha0=1.3;
estAlpha = fmincon(@logpdf,[estAlpha0],[],[],[],[],1,[],[],specoption);
totLL=-logpdf(estAlpha);
%     function [y] = MLEalpha

%     end

%calculate the log-likelihood
% if estAlpha==1
%     logconst=-log(log(xf)-log(xi));
% else
%     pp=1-estAlpha;
%     constant=-pp/(xi^pp-xf^pp);
%     logconst=log(-pp)-log(xi^pp-xf^pp);
% end
% lnpx=logconst-estAlpha*log(distributdata);


%calculate the Akaike information criterion, with small sample adjustment
n_params=1;   %numel(params);
LL=-h.*log(Prfunc(estAlpha,l,u));
AIC=-2*totLL+2*n_params*numel(distributdata)/(numel(distributdata)-n_params-1);

%%calculate quantities to carry out G-test
[Obs Expect] = GetObsExpect_AL(distributdata,'PowerLaw_bound',estAlpha,1);
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test
% plot(Obs,'+'); hold on plot(Expect); hold off
if plotting == 1
    fprintf('alpha-parameter : %5.4f \n',estAlpha)
    fprintf('totLL = %10.3f \n',totLL)
    fprintf('AIC = %10.3f \n',AIC)
    fprintf('pval = %5.4f \n',pval)
    figure; loglog((l+u)./2,h,'o')
    hold on
    loglog((l2+u2)/2,sum(h)*(Prfunc(estAlpha,l2,u2)),'r','LineWidth',2)
end
%% log pdf function
    function yfunc=logpdf(A)
        Pr=-log(xi^(1-A)-xf^(1-A))+log(l.^(1-A)-u.^(1-A));
        yfunc=-sum(h.*Pr);
    end
    function Prvals=Prfunc(A,lb,ub)
        Pr=-log(xi^(1-A)-xf^(1-A))+log(lb.^(1-A)-ub.^(1-A));
        Prvals=exp(Pr);
    end
end

