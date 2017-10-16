function [plambda AIC gof pval LL] = StrExpo_5Expo(distributdata,plotting,binway)
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

xmin=min(distributdata); xmax=Inf; % xmax=max(log(distributdata));
% distributdata(distributdata<xmin|distributdata>xmax)=[];   %exclude data that are outside bounds

U=unique(distributdata);
interval=U(2)-U(1);

%%set variables
k=6; %%order of exponential
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

g=@(a,a1,a2,a3,a4,a5,a6,l1,l2,l3,l4,l5,l6,x)(a1*exp(-l1*(exp(x)).^a)+a2*l2*exp(-l2*exp(x))+a3*l3*exp(-l3*exp(x))+a4*l4*exp(-l4*exp(x))+a5*l5*exp(-l5*exp(x))+a6*l6*exp(-l6*exp(x))).*exp(x);

miny=-1; minlambda=1./1500;
maxy=1; maxlambda=10;

% miny=-2; minlambda=0.99/max(x);
% maxy=0; maxlambda=0.99/min(x);
startpoint=[];
c=0;
while isempty(startpoint) && c<=2
    Allstartvals1=loghistfit_test(distributdata,0.2,6,300,0,0,1);
    Allstartvals2=loghistfit_test(distributdata,0.2,6,300,0,0,2);
    Allstartvals=[Allstartvals1;Allstartvals2];
    startpoint=Allstartvals(sum(Allstartvals(:,2:k),2)<1,:);
    c=c+1;
end
startpoint(:,k+1:end)=log(startpoint(:,k+1:end));
% Allstartvals(:,2:k)=normrnd(0,1,60,k-1);
% Allstartvals= RandInit(Numtrials,k,[miny,minlambda],[maxy,0.5*maxlambda]);
% Allstartvals(:,1)=[]; %%get rid of 1 set of y vals


specoption=optimset('Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',10000,'TolX',10^-6,'TolFun',10^-8,'Display','off','FinDiffType','central');
lowerbound=[0.1 miny miny miny miny miny log(10^-6) log(minlambda) log(minlambda) log(minlambda) log(minlambda) log(minlambda)];
upperbound=[1 maxy maxy maxy maxy maxy log(10) log(maxlambda) log(maxlambda) log(maxlambda) log(maxlambda) log(maxlambda)];
oparm=NaN(70,2*k); totLL=NaN(70,1);

parfor gg=1:min([70,size(startpoint,1)])
    
    try
        oparm(gg,:)=fmincon(@(A) str_5expo_pdf(distributdata,bin_bound,A),startpoint(gg,:),[0,-1,-1,-1,-1,-1,0,0,0,0,0,0;0,1,1,1,1,1,0,0,0,0,0,0],[0;2],[],[],lowerbound,upperbound,[],specoption);
        totLL(gg,1)=-str_5expo_pdf(distributdata,bin_bound,oparm(gg,:));
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
LL=-h.*log(Prfunc(oparm,l,u));
AIC=-2*totLL+2*n_params*numel(distributdata)/(numel(distributdata)-n_params-1);
%fprintf('AIC2 is %5.5f\n',AIC)

gof=[]; pval=[];
%%calculate quantities to carry out G-test
[Obs Expect] = GetObsExpect_AL(distributdata,'Str5Exponential',plambda,1);
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test %do the G-test

if plotting == 1
    fprintf('alpha-parameter : %5.4f \n',plambda(1))
    fprintf('a-parameters : %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n',plambda(2:k+1))
    fprintf('lambda-parameters : %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \t %5.4f \n',plambda(k+2:end))
    fprintf('totLL = %10.3f \n',totLL)
    fprintf('AIC = %10.3f \n',AIC)
    fprintf('pval = %5.4f \n',pval)
    figure; loglog((l+u)./2,h,'o')
    hold on 
    loglog((l2+u2)/2,sum(h)*(Prfunc(oparm(1,:),l2,u2)),'r','LineWidth',2)

%     loglog((uniques),(numUnique)/(sum(numUnique.*[diff(uniques); 0])),'ro',(uniques),...
%         exp(pdfit),'b')
    title('Str expo and 5 expo')
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
%     plot((edges(1:end-1)+edges(2:end))./2,g(plambda(1),histplot_A(1),histplot_A(2),histplot_A(3),histplot_A(4),histplot_A(5),histplot_A(6),plambda(8),plambda(9),plambda(10),plambda(11),plambda(12),plambda(13),(edges(1:end-1)+edges(2:end))./2).*diff(edges),'b')
%     hold off
end
%%%=======================================================
%%calculate the log(pdf) and return to 'mle' above
%     function fvals=str_5expo_pdf(A)
% 
%         alpha=A(1); a1=A(2); a2=A(3); a3=A(4); a4=A(5); a5=A(6); b=exp(A(7)); lam1=exp(A(8)); lam2=exp(A(9)); lam3=exp(A(10)); lam4=exp(A(11)); lam5=exp(A(12));
%         %calculate expo term
%         stexpo_term=log(a1)+log(gamma_incomplete(b*l.^alpha,1/alpha)-gamma_incomplete(b*u.^alpha,1/alpha))-log(gamma_incomplete(b*xmin^alpha,1/alpha));
%         expo_term1=log(a2)+lam1*xmin+log(exp(-lam1*l)-exp(-lam1*u));
%         expo_term2=log(a3)+lam2*xmin+log(exp(-lam2*l)-exp(-lam2*u));
%         expo_term3=log(a4)+lam3*xmin+log(exp(-lam3*l)-exp(-lam3*u));
%         expo_term4=log(a5)+lam4*xmin+log(exp(-lam4*l)-exp(-lam4*u));
%         expo_term5=log(1-a1-a2-a3-a4-a5)+lam5*xmin+log(exp(-lam5*l)-exp(-lam5*u));
%         stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
%         expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
%         expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
%         expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
%         expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
%         expo_term5 = reshape(expo_term5, 1, numel(expo_term5));
%         Pr=log(exp(expo_term5)+exp(expo_term4)+exp(expo_term3)+exp(expo_term2)+exp(expo_term1)+exp(stexpo_term));
%         fvals=-sum(h.*Pr);
%     end
    function PrVals=Prfunc(A,lb,ub)
        
        alpha=A(1); a1=A(2); a2=A(3); a3=A(4); a4=A(5); a5=A(6); b=exp(A(7)); lam1=exp(A(8)); lam2=exp(A(9)); lam3=exp(A(10)); lam4=exp(A(11)); lam5=exp(A(12));
        %calculate expo term
        stexpo_term=a1*(gamma_incomplete(b*lb.^alpha,1/alpha)-gamma_incomplete(b*ub.^alpha,1/alpha))./(gamma_incomplete(b*xmin^alpha,1/alpha));
        expo_term1=a2*exp(lam1*xmin)*(exp(-lam1*lb)-exp(-lam1*ub));
        expo_term2=a3*exp(lam2*xmin)*(exp(-lam2*lb)-exp(-lam2*ub));
        expo_term3=a4*exp(lam3*xmin)*(exp(-lam3*lb)-exp(-lam3*ub));
        expo_term4=a5*exp(lam4*xmin)*(exp(-lam4*lb)-exp(-lam4*ub));
        expo_term5=(1-a1-a2-a3-a4-a5)*exp(lam5*xmin)*(exp(-lam5*lb)-exp(-lam5*ub));
        stexpo_term = reshape(stexpo_term, 1, numel(stexpo_term));
        expo_term1 = reshape(expo_term1, 1, numel(expo_term1));
        expo_term2 = reshape(expo_term2, 1, numel(expo_term2));
        expo_term3 = reshape(expo_term3, 1, numel(expo_term3));
        expo_term4 = reshape(expo_term4, 1, numel(expo_term4));
        expo_term5 = reshape(expo_term5, 1, numel(expo_term5));
        PrVals=((expo_term5)+(expo_term4)+(expo_term3)+(expo_term2)+(expo_term1)+(stexpo_term));
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


