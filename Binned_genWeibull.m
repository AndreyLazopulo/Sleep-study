function [alph bet lamd AIC pval LL] = Binned_genWeibull(distributdata,plotting,binway)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    binway=1;
elseif nargin<2
    plotting=0;
    binway=1;
end
alph=[]; bet=[]; lamd=[]; AIC=[];
w = warning ('off','all');
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

meanx=mean(x); meanlnx=mean(log(x));  %%quantities that will be used later
specoption=optimset('Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',10000,'TolX',10^-6,'TolFun',10^-8,'Display','off','FinDiffType','central');

%%%initial values of alpha,a.k.a 'a' and beta, a.k.a 'b' come from the
%%%approximate estimations of Gamma distribution on Wikipedia
s=log(meanx)-meanlnx;
a0=0.6;%(3-s+sqrt(24*s+(s-3)^2))/(12*s);
b0=2.5;%meanx/a0;
lam0=0.4;%1;
% [param nLL]= fminsearch(@negloglike, [a0 b0 lam0]);
% param= mle(x,'logpdf',@logpdf,'logsf',@logsf,'start',[a0 b0 lam0],'optimfun','fmincon','options', otheroptions2);
param=fmincon(@logpdf,[a0 b0 lam0],[],[],[],[],[0 0 0],[1 Inf Inf],[],specoption);
alph=param(1); bet=param(2); lamd=param(3);
% totLL=-numel(x)*nLL;
totLL=logpdf(param);

%calculate the Akaike information criterion, with small sample adjustment
n_params=3;
LL=-h.*log(Prfunc(param,l,u));
AIC=2*totLL+2*n_params*numel(x)/(numel(x)-n_params-1);

%%calculate quantities to carry out G-test
[Obs Expect] = GetObsExpect_AL(x,'genWeibull',[alph,bet,lamd],1);
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test

% [Obs2 Expect2] = GetObsExpect_AL(x,'genWeibull',[0.2974,0.5316,0.4664],1);
% [gof pval2] = Gtest(Obs2,Expect2,n_params); %do the G-test

if plotting == 1
    fprintf('alpha parameter : %5.4f\n',alph)
    fprintf('betta parameter : %5.4f\n',bet)
    fprintf('lambda parameter : %5.4f\n',lamd)
    fprintf('totLL = %10.3f \n',totLL)
%     fprintf('totLL2 = %10.3f \n',logpdf([0.2974,0.5316,0.4664]))
    fprintf('AIC = %10.3f \n',AIC)
    fprintf('pval = %5.4f \n',pval)
%     fprintf('pval2 = %5.4f \n',pval2)
    figure; loglog((l+u)./2,h,'o')
    hold on 
    loglog((l2+u2)/2,sum(h)*(Prfunc(param,l2,u2)),'r','LineWidth',2)
    title('Generalized Weibull')
end
% figure;
% [uniques,numUnique] = count_unique(x);
% pdfit=logpdf(uniques,alph,bet,lamd);
% plot(log(uniques),log(numUnique),'ro',log(uniques),...
%     log(sum(numUnique.*[diff(uniques); 0]))+pdfit,'b')
%%%%============================================
    %% log pdf function
    function [mylogpdfvals]=logpdf(A)
        a=A(1); b=A(2); lam=A(3);
%         nn=(a+lam-1)/lam; 
%         zmin=double(vpa(expint(sym(nn),(xmin/b)^lam),5));
%         zmax=double(vpa(expint(sym(nn),(xmax/b)^lam),5));
%         zl=double(vpa(expint(sym(nn),(l./b).^lam),5));
%         zu=double(vpa(expint(sym(nn),(u./b).^lam),5));
        Pr=log(gamma_incomplete((l./b).^lam,(1-a)/lam)-gamma_incomplete((u./b).^lam,(1-a)/lam))-log(gamma_incomplete((xmin/b).^lam,(1-a)/lam)-gamma_incomplete((xmax/b).^lam,(1-a)/lam));
%         Pr=log(l.^(1-a).*zl-u.^(1-a).*zu)-log(xmin^(1-a)*zmin-xmax^(1-a)*zmax);
        mylogpdfvals=-sum(h.*Pr);
%         exint=xmin^(1-a)*zmin-xmax^(1-a)*zmax;
%         mylogpdfvals=-log(exint)-a*log(x)-(x./b).^lam+log(lam);
    end
    function [Prvals]=Prfunc(A,lb,ub)
        a=A(1); b=A(2); lam=A(3);
%         nn=(a+lam-1)/lam; 
%         zmin=double(vpa(expint(sym(nn),(xmin/b)^lam),5));
%         zmax=double(vpa(expint(sym(nn),(xmax/b)^lam),5));
%         zl=double(vpa(expint(sym(nn),(lb/b).^lam),25));
%         zr=double(vpa(expint(sym(nn),(ub/b).^lam),25));
        Pr=log(gamma_incomplete((lb./b).^lam,(1-a)/lam)-gamma_incomplete((ub./b).^lam,(1-a)/lam))-log(gamma_incomplete((xmin/b).^lam,(1-a)/lam)-gamma_incomplete((xmax/b).^lam,(1-a)/lam));
        Prvals=exp(Pr);
%         exint=xmin^(1-a)*zmin-xmax^(1-a)*zmax;
%         mylogpdfvals=-log(exint)-a*log(x)-(x./b).^lam+log(lam);
    end

end

