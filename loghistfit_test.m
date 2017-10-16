function [output,errors]=loghistfit_test(data,step,n,NtMax,zerorem,plotting,bintype)

% [filename,pathname]=uigetfile('.txt', 'Select test results files','MultiSelect','off');
% file1=fullfile(pathname,filename);
X=data; xmin=min(data); xmax=max(data);
% X=dlmread(file1);

[X_hist,Y_hist]=histcounts(X,min(X):1/3:600,'Normalization','probability');
[linfit,~]=fit(log((Y_hist(2:11)'+Y_hist(1:10)')./2),log(-log(X_hist(1:10)')),'poly1');
% figure; plot(linfit,log((Y_hist(2:end)'+Y_hist(1:end-1)')./2),log(-log(X_hist')))
fitpar=coeffvalues(linfit);
alpha_guess=fitpar(1);

N=length(X);
alimit=1;

if nargin < 7
    bintype = 1;
elseif nargin < 6
    bintype = 1;
    plotting = 0;
elseif nargin < 5
    bintype = 1;
    plotting = 0;
    zerorem = 0;
end

switch bintype
    case 1
        binbound=[log(unique(X));[max(log(unique(X)))+0.1:0.1:10]'];
    case 2
        binbound=[0:0.2:10]';
end

[fitY,edges]=histcounts(log(X),binbound);
fitY=fitY./length(X);
fitX=(edges(1:end-1)+edges(2:end))./2;
% fitY(5:2:9)=[];
% fitX(5:2:9)=[];
% fitY=h.Values;
m=zerorem;
I=find(fitY==0);
fitX(I(1:m))=[];
fitY(I(1:m))=[];
if plotting == 1
    figure;
    hold on
    histogram(log(X),edges,'Normalization','probability');
    xlim([0,10])
    plot(fitX,fitY,'o','Color','g')
end
fitY=fitY'./diff(edges);
% c=1;
% while c~=27
%     [Xdel,~,c]=ginput(1);
%     if c==1
%         [~,Idel]=min(abs(fitX-Xdel));
%         plot(fitX(Idel),fitY(Idel),'x','Color','r')
%         fitX(Idel)=[];
%         fitY(Idel)=[];
%     end
% end
%%Generate initial values for the parameters
% alpha0hat=1-N/(N*log(min(X))-sum(log(X)));
% alpha_Init=exprnd(alpha0hat/2,[NtMax,1]);  %%%random guesses for Alpha
% alpha_Init(alpha_Init>3)=3;  %%truncate if necessary

a_Init=normrnd(0,0.5,NtMax,n-1);  %%%random guesses for a's
a_Init(a_Init<-alimit)=-alimit; %%%%truncate values
a_Init(a_Init>alimit)=alimit;

minlam=10^-6; maxlam=1;
lamda_Init=exp(log(minlam)+(log(maxlam)-log(minlam))*rand(NtMax,n-1));
lamda_Init(lamda_Init>maxlam)=maxlam;
lamda_Init(lamda_Init<minlam)=minlam;
Allstartvals=[alpha_guess*2.5*ones(NtMax,1) a_Init 0.6*ones(NtMax,1) 0.1496*ones(NtMax,n-1)];

switch n
    case 1
        g=@(a,l1,x) exp(x).*exp(-l1*(exp(x)).^a)./(l1^(-1/a)*(gamma_incomplete(l1*xmin^a,1/a)-gamma_incomplete(l1*xmax^a,1/a))/a);
    case 2
        g=@(a,a1,l1,l2,x) exp(x).*(a1*exp(-l1*(exp(x)).^a)./(l1^(-1/a)*(gamma_incomplete(l1*xmin^a,1/a)-gamma_incomplete(l1*xmax^a,1/a))/a)+...
            (1-a1)*l2*exp(-l2*exp(x))./(exp(-l2*xmin)-exp(-l2*xmax)));
    case 3
        g=@(a,a1,a2,l1,l2,l3,x)(a1*exp(-l1*(exp(a*x)))./(l1^(-1/a)*(gamma_incomplete(l1*xmin^a,1/a)-gamma_incomplete(l1*xmax^a,1/a))/a)+a2*l2*exp(-l2*exp(x))./(exp(-l2*xmin)-exp(-l2*xmax))+(1-a1-a2)*l3*exp(-l3*exp(x))./(exp(-l3*xmin)-exp(-l3*xmax))).*exp(x);
    case 4
        g=@(a,a1,a2,a3,l1,l2,l3,l4,x)(a1*exp(-l1*(exp(a*x)))./(l1^(-1/a)*(gamma_incomplete(l1*xmin^a,1/a)-gamma_incomplete(l1*xmax^a,1/a))/a)+a2*l2*exp(-l2*exp(x))./(exp(-l2*xmin)-exp(-l2*xmax))+a3*l3*exp(-l3*exp(x))./(exp(-l3*xmin)-exp(-l3*xmax))+(1-a1-a2-a3)*l3*exp(-l4*exp(x))./(exp(-l4*xmin)-exp(-l4*xmax))).*exp(x);
    case 5
        g=@(a,a1,a2,a3,a4,l1,l2,l3,l4,l5,x)(a1*exp(-l1*(exp(a*x)))./(l1^(-1/a)*(gamma_incomplete(l1*xmin^a,1/a)-gamma_incomplete(l1*xmax^a,1/a))/a)+a2*l2*exp(-l2*exp(x))./(exp(-l2*xmin)-exp(-l2*xmax))+a3*l3*exp(-l3*exp(x))./(exp(-l3*xmin)-exp(-l3*xmax))+a4*l4*exp(-l4*exp(x))./(exp(-l4*xmin)-exp(-l4*xmax))+(1-a1-a2-a3-a4)*l5*exp(-l5*exp(x))./(exp(-l5*xmin)-exp(-l5*xmax))).*exp(x);
    case 6
        g=@(a,a1,a2,a3,a4,a5,l1,l2,l3,l4,l5,l6,x)(a1*exp(-l1*(exp(a*x)))./(l1^(-1/a)*(gamma_incomplete(l1*xmin^a,1/a)-gamma_incomplete(l1*xmax^a,1/a))/a)+a2*l2*exp(-l2*exp(x))./(exp(-l2*xmin)-exp(-l2*xmax))+a3*l3*exp(-l3*exp(x))./(exp(-l3*xmin)-exp(-l3*xmax))+a4*l4*exp(-l4*exp(x))./(exp(-l4*xmin)-exp(-l4*xmax))+a5*l5*exp(-l5*exp(x))./(exp(-l5*xmin)-exp(-l5*xmax))+(1-a1-a2-a3-a4-a5)*l6*exp(-l6*exp(x))./(exp(-l6*xmin)-exp(-l6*xmax))).*exp(x);
end
% astart=[5,-10,6,-8,4,-1];
% lstart=[0.5,0.01,0.05,0.1,0.07,0.3];
lb=[alpha_guess, -0.9*ones(1,n-1), 0.1, 0*ones(1,n-1)];
ub=[max([4,alpha_guess*4]), 0.9*ones(1,n-1), 10, 1*ones(1,n-1)];
% start=[2,astart(1:n),lstart(1:n)];
% preparms=zeros(NtMax,2*n);
parms=zeros(NtMax,2*n);
errors=zeros(NtMax,1);
r2=zeros(NtMax,1);
parfor gg=1:NtMax
    A0=Allstartvals(gg,:);
    try 
        [func,err]=fit(fitX,fitY,g,'Lower',lb,'Upper',ub,'Start',A0);
    catch
        continue
    end
    parms(gg,:)=coeffvalues(func);
    %     [xx indx]=sort(preparms(gg,n+1:end),'ascend');
    %     parms(gg,n+1:2*n)=xx;  %%%sorted exponents
    %     parms(gg,1:n)=preparms(gg,indx);   %%%sorted a's
    errors(gg)=err.rmse;
%     r2(gg)=err.adjrsquare;
end
if n==5
    for kk=1:length(parms)
        oparm=[parms(kk,1:n),log(parms(kk,n+1:end))];
        if ~isreal(-str_4expo_pdf(X,[unique(X);Inf],oparm));
            parms(kk,:)=NaN(1,2*n);
            errors(kk)=NaN;
        end
    end
elseif n==6
    for kk=1:length(parms)
        oparm=[parms(kk,1:n),log(parms(kk,n+1:end))];
        if ~isreal(-str_5expo_pdf(X,[unique(X);Inf],oparm));
            parms(kk,:)=NaN(1,2*n);
            errors(kk)=NaN;
        end
    end
end
% parms=[parms(:,1:n-1),1-sum(parms(:,1:n-1),2),parms(:,n:2*n-1)];
[errors,outind]=sort(errors);
rounderr=round(10^4*errors)/10^4;
% ind=find(rounderr==min(rounderr));
parms=parms(outind,:);
parms=parms(1:100,:);
errors=rounderr(1:100,:);
% parms=parms(rounderr==min(rounderr),:);
if plotting == 1
    hold on
    switch n
        case 1
            plot(fitX,g(parms(1,1),parms(1,2),fitX).*diff(edges'),'r')
        case 2
            plot(fitX,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),fitX).*diff(edges),'r')
        case 3
            plot(fitX,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),fitX).*diff(edges),'r')
        case 4
            plot(fitX,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),parms(1,7),parms(1,8),fitX).*diff(edges),'r')
        case 5
            plot(fitX,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),parms(1,7),parms(1,8),parms(1,9),parms(1,10),fitX).*diff(edges),'r')
        case 6
            plot(fitX,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),parms(1,7),parms(1,8),parms(1,9),parms(1,10),parms(1,11),parms(1,12),fitX).*diff(edges),'r')
    end
end
output=[parms(1:40,1:n),parms(1:40,n+1:end)];
% for gg=1:60
%     output(gg,2:n)=parms(gg,2:n)./sum(parms(gg,2:n+1));
%     output(gg,2:n)=log(output(gg,2:n)./(1-sum(output(gg,2:n))));
% end
% for v=1:100
%     for j=1:n-1
%             lamda=(parms(v,j+n));
%             norm=exp(-lamda*xmin)-exp(-lamda*xmax);
%             parms(v,j)=parms(v,j)*norm;
%     end
% end
% output=zeros(50,2*n-1);
% j=1;
% jj=1;
% while jj<=50
%     const=1./(1-sum(parms(j,1:n-1)));
%     if const>0
%         output(jj,:)=[log(parms(j,1:n-1)*const),parms(j,n+1:2*n)];
%         jj=jj+1;
%         j=j+1;
%     else
%         j=j+1;
%     end
% end
end
