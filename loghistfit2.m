function [parms,errors,r2]=loghistfit2(data,step,n,NtMax,zerorem,plotting)
% [filename,pathname]=uigetfile('.txt', 'Select test results files','MultiSelect','off');
% file1=fullfile(pathname,filename);
X=data;
% X=dlmread(file1);
N=length(X);
alimit=1;
[fitY,~]=histcounts(log(X),0:step:10);
fitY=fitY./length(X);
fitX=step/2:step:10;
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
    histogram(log(X),-10:step:10,'Normalization','probability');
    xlim([0,10])
    plot(fitX,fitY,'o','Color','g')
end
fitY=fitY./step;
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

a_Init=randi([0,alimit],NtMax,n);  %%%random guesses for a's
a_Init(a_Init<-alimit)=0; %%%%truncate values
a_Init(a_Init>alimit)=alimit;

minlam=10^-6; maxlam=1;
lamda_Init=exp(log(minlam)+(log(maxlam)-log(minlam))*rand(NtMax,n));
lamda_Init(lamda_Init>maxlam)=maxlam;
lamda_Init(lamda_Init<minlam)=minlam;
Allstartvals=[a_Init lamda_Init];

switch n
    case 2
        g=@(a1,a2,l1,l2,x)(a1*l1*exp(-l1*exp(x))+a2*l2*exp(-l2*exp(x))).*exp(x);
    case 3
        g=@(a1,a2,a3,l1,l2,l3,x)(a1*l1*exp(-l1*exp(x))+a2*l2*exp(-l2*exp(x))+a3*l3*exp(-l3*exp(x))).*exp(x);
    case 4
        g=@(a1,a2,a3,a4,l1,l2,l3,l4,x)(a1*l1*exp(-l1*exp(x))+a2*l2*exp(-l2*exp(x))+a3*l3*exp(-l3*exp(x))+a4*l4*exp(-l4*exp(x))).*exp(x);
    case 5
        g=@(a1,a2,a3,a4,a5,l1,l2,l3,l4,l5,x)(a1*l1*exp(-l1*exp(x))+a2*l2*exp(-l2*exp(x))+a3*l3*exp(-l3*exp(x))...
            +a4*l4*exp(-l4*exp(x))+a5*l5*exp(-l5*exp(x))).*exp(x);
    case 6
        g=@(a1,a2,a3,a4,a5,a6,l1,l2,l3,l4,l5,l6,x)(a1*l1*exp(-l1*exp(x))+a2*l2*exp(-l2*exp(x))+a3*l3*exp(-l3*exp(x))...
            +a4*l4*exp(-l4*exp(x))+a5*l5*exp(-l5*exp(x))+a6*l6*exp(-l6*exp(x))).*exp(x);
end
% astart=[5,-10,6,-8,4,-1];
% lstart=[0.5,0.01,0.05,0.1,0.07,0.3];
lb=[repmat(0,1,n),zeros(1,n)];
ub=[repmat(1,1,n),ones(1,n)];
% start=[2,astart(1:n),lstart(1:n)];
% preparms=zeros(NtMax,2*n);
parms=zeros(NtMax,2*n);
errors=zeros(NtMax,1);
r2=zeros(NtMax,1);
for gg=1:NtMax
    A0=Allstartvals(gg,:);
    [func,err]=fit(fitX',fitY',g,'Lower',lb,'Upper',ub,'Start',A0);
    parms(gg,1:2*n)=coeffvalues(func);
    %     [xx indx]=sort(preparms(gg,n+1:end),'ascend');
    %     parms(gg,n+1:2*n)=xx;  %%%sorted exponents
    %     parms(gg,1:n)=preparms(gg,indx);   %%%sorted a's
    errors(gg)=err.rmse;
    r2(gg)=err.adjrsquare;
end
% parms=[parms(:,1:n-1),1-sum(parms(:,1:n-1),2),parms(:,n:2*n-1)];
[errors,outind]=sort(errors);
rounderr=round(10^4*errors)/10^4;
% ind=find(rounderr==min(rounderr));
parms=parms(outind,:);
parms=parms(1:30,:);
errors=rounderr(1:30,:);
% parms=parms(rounderr==min(rounderr),:);
if plotting == 1
    hold on
    switch n
        case 2
            plot(0:0.1:10,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),0:0.1:10)*step,'r')
        case 3
            plot(0:0.1:10,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),0:0.1:10)*step,'r')
        case 4
            plot(0:0.1:10,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),parms(1,7),parms(1,8),0:0.1:10)*step,'r')
        case 5
            plot(0:0.1:10,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),parms(1,7),parms(1,8),parms(1,9),parms(1,10),0:0.1:10)*step,'r')
        case 6
            plot(0:0.1:10,g(parms(1,1),parms(1,2),parms(1,3),parms(1,4),parms(1,5),parms(1,6),parms(1,7),parms(1,8),parms(1,9),parms(1,10),parms(1,11),parms(1,12),0:0.1:10)*step,'r')
    end
end
end
