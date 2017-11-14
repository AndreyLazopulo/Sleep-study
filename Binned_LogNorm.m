function [mu sigma AIC gof pval LL]=Binned_LogNorm(distributdata,plotting,binway)
%%%Function computes the maximum log-likelihood values for a log-normal parameters to the input
%%%distribution Sheyum. 07/05/13
%%%modified Lazopulo 12/28/16
if nargin<3
    binway=1;
elseif nargin<2
    plotting=0;
    binway=1;
end
mu=[]; sigma=[]; AIC=[]; gof=[]; pval=[];
xmin=min(distributdata);
xmax=max(distributdata);

U=unique(distributdata);
interval=U(2)-U(1);

switch binway
    case 1
        bin_bound=[unique(distributdata);xmax+interval];
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

x=distributdata;
x(x<xmin|x>xmax)=[];   %exclude data outside range
specoption=optimset('Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',10000,'TolX',10^-6,'TolFun',10^-8,'Display','off','FinDiffType','central');
%%initial guesses for the two parameters
mu0=mean(log(x));
sig0=sqrt(sum((log(x)-mu0).^2)/numel(x));
%%%do minimization of the negative of the totLL kernel
% parm= mle(x,'logpdf',@logpdf,'start',[mu0 sig0],'optimfun','fmincon','options', otheroptions2);
parm=fmincon(@logpdf,[mu0,sig0],[],[],[],[],[],[],[],specoption);
totLL=-logpdf(parm);
mu=parm(1); sigma=parm(2);  %%max likelihood estimators

%calculate the total likelihood
% totLL=-numel(x)*(mean(log(x))+knLL-log(sqrt(2/pi)));

%calculate the Akaike information criterion, with small sample adjustment
n_params=2;
LL=-h.*log(Prfunc(parm,l,u));
AIC=-2*totLL+2*n_params*numel(x)/(numel(x)-n_params-1);

%%calculate quantities to carry out G-test
[Obs Expect] = GetObsExpect_AL(x,'Log-Norm',[mu,sigma],1);
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test

if plotting == 1
    fprintf('mu parameter : %5.4f\n',mu)
    fprintf('sigma parameter : %5.4f\n',sigma)
    fprintf('totLL = %10.3f \n',totLL)
    fprintf('AIC = %10.3f \n',AIC)
    fprintf('pval = %5.4f \n',pval)
    figure; loglog((l+u)./2,h,'o')
    hold on 
    loglog((l2+u2)/2,sum(h)*(Prfunc(parm,l2,u2)),'r','LineWidth',2)
end
% figure;
% [uniques,numUnique] = count_unique(x);
% pdfit=logpdf(uniques,mu,sigma);
% plot(log(uniques),log(numUnique),'ro',log(uniques),...
%     log(sum(numUnique.*[diff(uniques); 0]))+pdfit,'b')
%% log pdf function
    function yfunc=logpdf(A)
        m=A(1); sig=A(2);
        deno=sqrt(2)*sig;
        zmin=(m-log(xmin))/deno;
        zmax=(m-log(xmax))/deno;
        zl=(m-log(l))./deno;
        zr=(m-log(u))./deno;
        Pr=-log(erf(zmin)-erf(zmax))+log(erf(zl)-erf(zr));
        yfunc=-sum(h.*Pr);
    end
    function PrVals=Prfunc(A,lb,ub)
        m=A(1); sig=A(2);
        deno=sqrt(2)*sig;
        zmin=(m-log(xmin))/deno;
        zmax=(m-log(xmax))/deno;
        zl=(m-log(lb))./deno;
        zr=((m-log(ub))./deno);
        Pr=-log(erf(zmin)-erf(zmax))+log(erf(zl)-erf(zr));
        PrVals=exp(Pr);
    end
    
end