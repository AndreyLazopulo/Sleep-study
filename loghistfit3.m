function [output,errors]=loghistfit3(data,step,n,NtMax,zerorem,plotting,bintype)
% [filename,pathname]=uigetfile('.txt', 'Select test results files','MultiSelect','off');
% file1=fullfile(pathname,filename);
X=data; xmin=min(data); xmax=max(data);
N=length(X);
% alimit=5;

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
%initial guess for alpha
% [Y,E]=histcounts(X,min(X):1/3:max(X)+1/3);
% f=fit(log(E(1:4)'),log(Y(1:4)'),'poly1');
% CC=coeffvalues(f);
alpha0=-5;%CC(1);

[fitY,edges]=histcounts(log(X),binbound);
fitY=fitY./length(X);
fitX=(edges(1:end-1)+edges(2:end))/2;
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
    xlim([0;10])
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
alpha_Init=0.4*alpha0*rand(NtMax,1)+alpha0-0.2*alpha0;

a_Init=normrnd(0,1,NtMax,n-1);  %%%random guesses for a's
a_Init(a_Init<-1)=-1; %%%%truncate values
a_Init(a_Init>1)=1;

minlam=1/1500; maxlam=10;
lamda_Init=exp(log(minlam)+(log(maxlam)-log(minlam))*rand(NtMax,n));
lamda_Init(lamda_Init>maxlam)=maxlam;
lamda_Init(lamda_Init<minlam)=minlam;
Allstartvals=[alpha_Init a_Init lamda_Init];

switch n
    case 1
        g=@(a,l1,x)(l1^(1-a)*exp(-l1*exp(x))./(-gamma_incomplete(xmax*l1,1-a)+gamma_incomplete(xmin*l1,1-a))).*exp((1-a)*x);
    case 2
        g=@(a,a1,l1,l2,x)(a1*l1^(1-a)*exp(-l1*exp(x))./(-gamma_incomplete(xmax*l1,1-a)+gamma_incomplete(xmin*l1,1-a))+(1-a1)*l2^(1-a)*exp(-l2*exp(x))./(-gamma_incomplete(xmax*l2,1-a)+gamma_incomplete(xmin*l2,1-a))).*exp((1-a)*x);
    case 3
        g=@(a,a1,a2,l1,l2,l3,x)(a1*l1^(1-a)*exp(-l1*exp(x))./(-gamma_incomplete(xmax*l1,1-a)+gamma_incomplete(xmin*l1,1-a))+a2*l2^(1-a)*exp(-l2*exp(x))./(-gamma_incomplete(xmax*l2,1-a)+gamma_incomplete(xmin*l2,1-a))...
            +(1-a1-a2)*l3^(1-a)*exp(-l3*exp(x))./(-gamma_incomplete(xmax*l3,1-a)+gamma_incomplete(xmin*l3,1-a))).*exp((1-a)*x);
    case 4
        g=@(a,a1,a2,a3,l1,l2,l3,l4,x)(a1*l1^(1-a)*exp(-l1*exp(x))./(-gamma_incomplete(xmax*l1,1-a)+gamma_incomplete(xmin*l1,1-a))+a2*l2^(1-a)*exp(-l2*exp(x))./(-gamma_incomplete(xmax*l2,1-a)+gamma_incomplete(xmin*l2,1-a))...
            +a3*l3^(1-a)*exp(-l3*exp(x))./(-gamma_incomplete(xmax*l3,1-a)+gamma_incomplete(xmin*l3,1-a))+(1-a1-a2-a3)*l4^(1-a)*exp(-l4*exp(x))./(-gamma_incomplete(xmax*l4,1-a)+gamma_incomplete(xmin*l4,1-a))).*exp((1-a)*x);
    case 5
        g=@(a,a1,a2,a3,a4,l1,l2,l3,l4,l5,x)(a1*l1^(1-a)*exp(-l1*exp(x))./(-gamma_incomplete(xmax*l1,1-a)+gamma_incomplete(xmin*l1,1-a))+a2*l2^(1-a)*exp(-l2*exp(x))./(-gamma_incomplete(xmax*l2,1-a)+gamma_incomplete(xmin*l2,1-a))...
            +a3*l3^(1-a)*exp(-l3*exp(x))./(-gamma_incomplete(xmax*l3,1-a)+gamma_incomplete(xmin*l3,1-a))+a4*l4^(1-a)*exp(-l4*exp(x))./(-gamma_incomplete(xmax*l4,1-a)+gamma_incomplete(xmin*l4,1-a))...
            +(1-a1-a2-a3-a4)*l5^(1-a)*exp(-l5*exp(x))./(-gamma_incomplete(xmax*l5,1-a)+gamma_incomplete(xmin*l5,1-a))).*exp((1-a)*x);
    case 6
        g=@(a,a1,a2,a3,a4,a5,l1,l2,l3,l4,l5,l6,x)(a1*l1^(1-a)*exp(-l1*exp(x))./(-gamma_incomplete(xmax*l1,1-a)+gamma_incomplete(xmin*l1,1-a))+a2*l2^(1-a)*exp(-l2*exp(x))./(-gamma_incomplete(xmax*l2,1-a)+gamma_incomplete(xmin*l2,1-a))...
            +a3*l3^(1-a)*exp(-l3*exp(x))./(-gamma_incomplete(xmax*l3,1-a)+gamma_incomplete(xmin*l3,1-a))+a4*l4^(1-a)*exp(-l4*exp(x))./(-gamma_incomplete(xmax*l4,1-a)+gamma_incomplete(xmin*l4,1-a))...
            +a5*l5^(1-a)*exp(-l5*exp(x))./(-gamma_incomplete(xmax*l5,1-a)+gamma_incomplete(xmin*l5,1-a))+(1-a1-a2-a3-a4-a5)*l6^(1-a)*exp(-l6*exp(x))./(-gamma_incomplete(xmax*l6,1-a)+gamma_incomplete(xmin*l6,1-a))).*exp((1-a)*x);
end
% astart=[5,-10,6,-8,4,-1];
% lstart=[0.5,0.01,0.05,0.1,0.07,0.3];
lb=[alpha0-0.2*alpha0,-1*ones(1,n-1),repmat(minlam,1,n)];
ub=[alpha0-0.2*alpha0,ones(1,n-1),10*ones(1,n)];
% start=[2,astart(1:n),lstart(1:n)];
parms=zeros(NtMax,2*n);
errors=zeros(NtMax,1);
parfor gg=1:NtMax
    A0=Allstartvals(gg,:);
    [func,err]=fit(fitX,fitY,g,'Lower',lb,'Upper',ub,'Start',A0);
    parms(gg,:)=coeffvalues(func);
    errors(gg)=err.rmse;
end
if n==5
    for kk=1:length(parms)
        oparm=[parms(kk,1:n),log(parms(kk,n+1:end))];
        if ~isreal(-binned_pl5expo_pdf(X,[unique(X);Inf],oparm)) || sum(oparm(2:n))>1
            parms(kk,:)=NaN(1,2*n);
            errors(kk)=NaN;
        end
    end
elseif n==6
    for kk=1:length(parms)
        oparm=[parms(kk,1:n),log(parms(kk,n+1:end))];
        if ~isreal(-binned_pl6expo_pdf(X,[unique(X);Inf],oparm)) || sum(oparm(2:n))>1
            parms(kk,:)=NaN(1,2*n);
            errors(kk)=NaN;
        end
    end
end

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
            plot(fitX,g(parms(1,1),parms(1,2),fitX).*diff(edges),'r')
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
% for v=1:100
%     for j=1:n-1
%             lamda=(parms(v,j+n+1));
%             norm=-gamma_incomplete(xmax*(lamda),1-parms(v,1))+...
%                 gamma_incomplete(xmin*(lamda),1-parms(v,1));
%             parms(v,j+1)=parms(v,j+1)*norm;
%     end
% end
output=parms(1:40,:);
% j=1;
% jj=1;
% while jj<=50
%     const=1/(1-sum(parms(j,2:n)));
%     if const>0
%         output(jj,:)=[log(parms(j,2:n)*const),parms(j,1),parms(j,n+2:2*n+1)];
%         jj=jj+1;
%         j=j+1;
%     else
%         j=j+1;
%     end
% end
