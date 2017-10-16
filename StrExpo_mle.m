function [plambda AIC gof pval LL] = StrExpo_mle(distributdata,plotting,binway)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%added sorting of parameters, (in increasing order) on 12/25/2013
if nargin<3
    plotting=0;
end
Numtrials=10;% (nargin>3)*Nt+(nargin<=3)*10;
plambda=[]; lnpx=[]; AIC=[]; gof=[]; pval=[];  %initialize

xmin=min(distributdata); xmax=Inf;% xmax=max(log(distributdata));
% distributdata(distributdata<xmin|distributdata>xmax)=[];   %exclude data that are outside bounds

U=unique(distributdata);
interval=U(2)-U(1);

%%set variables
k=1; %%order of exponential
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

g=@(a,a1,l1,x) a1*exp(-l1*(exp(x)).^a).*exp(x);

miny=0.001; minlambda=1./1500;
maxy=1; maxlambda=1./1;

% miny=-2; minlambda=0.99/max(x);
% maxy=0; maxlambda=0.99/min(x);
Allstartvals1=loghistfit_test(distributdata,0.2,1,300,0,0,1);
Allstartvals2=loghistfit_test(distributdata,0.2,1,300,0,0,2);
Allstartvals=[Allstartvals1(1:30,:);Allstartvals2(1:30,:)];
Allstartvals(:,2)=log(Allstartvals(:,2));
% Allstartvals(:,2:k)=normrnd(0,1,60,k-1);
% Allstartvals= RandInit(Numtrials,k,[miny,minlambda],[maxy,0.5*maxlambda]);
% Allstartvals(:,1)=[]; %%get rid of 1 set of y vals


specoption=optimset('Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',10000,'TolX',10^-6,'TolFun',10^-8,'Display','off','FinDiffType','central');
lowerbound=[0.1 log(10^-6)];
upperbound=[1 log(10)];
oparm=NaN(50,2); totLL=NaN(50,1);

for gg=1:min(50,size(Allstartvals,1))
    try
        oparm(gg,:)=fmincon(@strexpo_pdf,Allstartvals(gg,:),[],[],[],[],lowerbound,upperbound,[],specoption);
        totLL(gg,1)=-strexpo_pdf(oparm(gg,:));
    catch
        continue
    end
    
end
mm=find(totLL==max(totLL),1,'first');
totLL=totLL(mm);
oparm=oparm(mm,:);
plambda=zeros(1,2*k+1); %%initialize
%%inverse transform of the prefactors in the exponential terms
plambda(1)=oparm(1);
plambda(2:k)=oparm(2:k); %exp(oparm(2:k))./(1+sum(exp(oparm(1:k))));
plambda(k+1)=1-sum(plambda(2:k));
%%save values of lambda
plambda(k+2:end)=exp(oparm(k+1:end));

%calculate the Akaike information criterion, with small sample adjustment
n_params=2;
LL=-h.*Prfunc(oparm,l,u);
AIC=-2*totLL+2*n_params*numel(distributdata)/(numel(distributdata)-n_params-1);
%fprintf('AIC2 is %5.5f\n',AIC)

gof=[]; pval=[];
%%calculate quantities to carry out G-test
[Obs Expect] = GetObsExpect_AL(distributdata,'StrExpo',[plambda(1),plambda(3)],1);
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test
% figure; plot(h,sum(h)*exp(Prfunc(oparm(1,:),l,u)),'o',h,h);
if plotting == 1
    fprintf('alpha-parameter : %5.4f \n',plambda(1))
    fprintf('a-parameters : %5.4f \n',plambda(2:k+1))
    fprintf('lambda-parameters : %5.4f \n',plambda(k+2:end))
    fprintf('totLL = %10.3f \n',totLL)
    fprintf('AIC = %10.3f \n',AIC)
    fprintf('pval = %5.4f \n',pval)
    figure; loglog((l+u)./2,h,'o')
    hold on 
    loglog((l2+u2)/2,sum(h)*exp(Prfunc(oparm(1,:),l2,u2)),'r','LineWidth',2)

%     loglog((uniques),(numUnique)/(sum(numUnique.*[diff(uniques); 0])),'ro',(uniques),...
%         exp(pdfit),'b')
    title('Str Expo')
%     edges=[log(unique(distributdata));[max(log(unique(distributdata)))+0.1:0.1:10]'];
%     b_par=plambda(2); a=(plambda(1));
%     norm2=b_par^(-1/a)*(gamma_incomplete(b_par*xmin^a,1/a)-gamma_incomplete(b_par*xmax^a,1/a))/a;
%     histplot_A(1)=1/norm2;
%     figure;
%     histogram(log(distributdata),edges,'Normalization','probability');
%     hold on
%     plot((edges(1:end-1)+edges(2:end))./2,g(plambda(1),histplot_A(1),plambda(2),(edges(1:end-1)+edges(2:end))./2).*diff(edges),'b')
%     hold off
end
%%%=======================================================
%%calculate the log(pdf) and return to 'mle' above
    function fvals=strexpo_pdf(A)

        alpha=A(1); b=exp(A(2)); a1=1;
        %calculate expo term
        stexpo_term=log(a1)+log(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))-log(gamma_incomplete(b*xmin^alpha,1/alpha));
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        Pr=stexpo_term;
        fvals=-sum(h.*Pr);
    end
    function PrVals=Prfunc(A,lb,ub)
        
        alpha=A(1); b=exp(A(2)); a1=1;
        %calculate expo term
        stexpo_term=log(a1)+log(gamma_incomplete(b*lb.^alpha,1/alpha)-gamma_incomplete(b*ub.^alpha,1/alpha))-log(gamma_incomplete(b*xmin^alpha,1/alpha));
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        PrVals=stexpo_term;
    end
%%calculate the log(sf) or log(1-cdf) and return to 'mle' above
%     function logsfvals=logsf(x,sy1,sl1,sl2)
%         pv=zeros(1,k-1);  terms=zeros(1); %%initialize variables
%         p=[sy1]; lb=[sl1 sl2];
%         %%transform p values
%         for j=1:(k-1),pv(j)=exp(p(j))/(1+sum(exp(p))); end
%         %calculate last term of distribution
%         terms=(1-sum(pv))*(exp(lb(k)*(xmax-x))-1);
%         terms=terms./(exp(lb(k)*(xmax-xmin))-1);
%         %%calculate the full distribution
%         for jj=1:(k-1)
%             constant=exp(lb(jj)*(xmax-xmin))-1;
%             newterm=pv(jj)*(exp(lb(jj)*(xmax-x))-1);
%             terms=terms+newterm./constant;
%         end
%         logsfvals=log(terms);
%     end

end


