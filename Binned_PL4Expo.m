function [plambda AIC gof pval LL] = Binned_PL4Expo(distributdata,plotting,binway)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%%added sorting of parameters, (in increasing order) on 12/25/2013
if nargin<3
    binway=1;
elseif nargin<2
    plotting=0;
    binway=1;
end
Numtrials=10;% (nargin>3)*Nt+(nargin<=3)*10;
plambda=[]; lnpx=[]; AIC=[]; gof=[]; pval=[];  %initialize

xmin=min(distributdata); xmax=Inf;% xmax=max(log(distributdata));
% distributdata(distributdata<xmin|distributdata>xmax)=[];   %exclude data that are outside bounds

U=unique(distributdata);
interval=U(2)-U(1);

%%set variables
k=4; %%order of exponential
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
% l2=1:1/3:max(distributdata); u2=l2+1/3;

g=@(a,a1,a2,a3,l1,l2,l3,x)(a1*l1^(1-a)*exp(-l1*exp(x))+a2*l2^(1-a)*exp(-l2*exp(x))+a3*l3^(1-a)*exp(-l3*exp(x))).*exp((1-a)*x);

miny=-1; minlambda=1./1500;
maxy=1; maxlambda=10;

% miny=-2; minlambda=0.99/max(x);
% maxy=0; maxlambda=0.99/min(x);
Allstartvals=[];
c=0;
while isempty(Allstartvals) && c<=2
    Allstartvals1=loghistfit3(distributdata,0.2,4,300,0,0,1);
    Allstartvals2=loghistfit3(distributdata,0.2,4,300,0,0,2);
    Allstartvals=[Allstartvals1(1:30,:);Allstartvals2(1:30,:)];
    Allstartvals=Allstartvals(sum(Allstartvals(:,2:k),2)<1,:);
    Allstartvals(:,k+1:end)=log(Allstartvals(:,k+1:end));
    c=c+1;
end

% Allstartvals(:,2:k)=normrnd(0,1,60,k-1);
% Allstartvals= RandInit(Numtrials,k,[miny,minlambda],[maxy,0.5*maxlambda]);
% Allstartvals(:,1)=[]; %%get rid of 1 set of y vals


specoption=optimset('Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',10000,'TolX',10^-6,'TolFun',10^-8,'Display','off','FinDiffType','central');
lowerbound=[0 miny miny miny log(minlambda) log(minlambda) log(minlambda) log(minlambda)];
upperbound=[3 maxy maxy maxy log(maxlambda) log(maxlambda) log(maxlambda) log(maxlambda)];
oparm=NaN(50,2*k); totLL=NaN(50,1);

parfor gg=1:min(50,size(Allstartvals,1))
    try
        oparm(gg,:)=fmincon(@(A) binned_pl4expo_pdf(distributdata,bin_bound,A),Allstartvals(gg,:),[0,-1,-1,-1,0,0,0,0;0,1,1,1,0,0,0,0],[0;2],[],[],lowerbound,upperbound,[],specoption);
        totLL(gg,1)=-binned_pl4expo_pdf(distributdata,bin_bound,oparm(gg,:));
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
n_params=2*k;
LL=-h.*Prfunc(oparm,l,u);
AIC=-2*totLL+2*n_params*numel(distributdata)/(numel(distributdata)-n_params-1);
%fprintf('AIC2 is %5.5f\n',AIC)

gof=[]; pval=[];
%%calculate quantities to carry out G-test
[Obs Expect] = GetObsExpect_AL(distributdata,'PL4Exponential',plambda,1);
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test
% figure; plot(h,sum(h)*exp(Prfunc(oparm(1,:),l,u)),'o',h,h);

if plotting == 1
    fprintf('alpha-parameter : %5.4f \n',plambda(1))
    fprintf('a-parameters : %5.4f \t %5.4f \t %5.4f \t %5.4f \n',plambda(2:k+1))
    fprintf('lambda-parameters : %5.4f \t %5.4f \t %5.4f \t %5.4f \n',plambda(k+2:end))
    fprintf('totLL = %10.3f \n',totLL)
    fprintf('AIC = %10.3f \n',AIC)
    fprintf('pval = %5.4f \n',pval)
    figure; loglog((l+u)./2,h,'o')
    hold on 
    loglog((l2+u2)/2,sum(h)*exp(Prfunc(oparm(1,:),l2,u2)),'r','LineWidth',2)

%     loglog((uniques),(numUnique)/(sum(numUnique.*[diff(uniques); 0])),'ro',(uniques),...
%         exp(pdfit),'b')
    title('PL with 4 Expo')
%     edges=[log(unique(distributdata));[max(log(unique(distributdata)))+0.1:0.1:10]'];
%     b_par=plambda(k+2); a=(plambda(1));
%     norm2=b_par^(-1/a)*(gamma_incomplete(b_par*xmin^a,1/a)-gamma_incomplete(b_par*xmax^a,1/a))/a;
%     histplot_A(1)=plambda(2)/norm2;
%     for jjj=2:k
%         lamda=plambda(k+1+jjj);
%         norm2=exp(-lamda*xmin)-exp(-lamda*xmax);
%         histplot_A(jjj)=plambda(jjj+1)/norm2;
%     end
%     figure;
%     histogram(log(distributdata),edges,'Normalization','probability');
%     hold on
%     plot((edges(1:end-1)+edges(2:end))./2,g(plambda(1),histplot_A(1),histplot_A(2),plambda(4),plambda(5),(edges(1:end-1)+edges(2:end))./2).*diff(edges),'b')
%     hold off
end
%%%=======================================================
%%calculate the log(pdf) and return to 'mle' above
%     function fvals=binned_pl4expo_pdf(A)
% 
%         alpha=A(1); a1=A(2); a2=A(3); a3=A(4); lam1=exp(A(5)); lam2=exp(A(6)); lam3=exp(A(7)); lam4=exp(A(8));
%         %calculate expo term
%         expo_term1=log(a1)+log(gamma_incomplete(lam1*l,1-alpha)-gamma_incomplete(lam1*u,1-alpha))-log(gamma_incomplete(lam1*xmin,1-alpha));
%         expo_term2=log(a2)+log(gamma_incomplete(lam2*l,1-alpha)-gamma_incomplete(lam2*u,1-alpha))-log(gamma_incomplete(lam2*xmin,1-alpha));
%         expo_term3=log(a3)+log(gamma_incomplete(lam3*l,1-alpha)-gamma_incomplete(lam3*u,1-alpha))-log(gamma_incomplete(lam3*xmin,1-alpha));
%         expo_term4=log(1-a1-a2-a3)+log(gamma_incomplete(lam4*l,1-alpha)-gamma_incomplete(lam4*u,1-alpha))-log(gamma_incomplete(lam4*xmin,1-alpha));
%         expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
%         expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
%         expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
%         expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
%         Pr=log(exp(expo_term1)+exp(expo_term2)+exp(expo_term3)+exp(expo_term4));
%         fvals=-sum(h.*Pr);
%     end
    function PrVals=Prfunc(A,lb,ub)
        
        alpha=A(1); a1=A(2); a2=A(3); a3=A(4); lam1=exp(A(5)); lam2=exp(A(6)); lam3=exp(A(7)); lam4=exp(A(8));
        %calculate expo term
        expo_term1=a1*(gamma_incomplete(lam1*lb,1-alpha)-gamma_incomplete(lam1*ub,1-alpha))./(gamma_incomplete(lam1*xmin,1-alpha));
        expo_term2=a2*(gamma_incomplete(lam2*lb,1-alpha)-gamma_incomplete(lam2*ub,1-alpha))./(gamma_incomplete(lam2*xmin,1-alpha));
        expo_term3=a3*(gamma_incomplete(lam3*lb,1-alpha)-gamma_incomplete(lam3*ub,1-alpha))./(gamma_incomplete(lam3*xmin,1-alpha));
        expo_term4=(1-a1-a2-a3)*(gamma_incomplete(lam4*lb,1-alpha)-gamma_incomplete(lam4*ub,1-alpha))./(gamma_incomplete(lam4*xmin,1-alpha));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        PrVals=log((expo_term1)+(expo_term2)+(expo_term3)+(expo_term4));
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


