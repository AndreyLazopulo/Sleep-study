function ccdfs = PlotDataModelCdf(empdistr, models, params,xmin,xmax,datName,figNum)

%function plots and returns the empirical and theoretical ccdfs 

ccdfs=[];
if nargin<7,
    figNum=1;
elseif nargin<6
    datname='Data'; figNum=1;
end

N_plots=size(models,1);
n1=3;
n2=ceil(N_plots/3);
cset=varycolor(N_plots);

%%%first plot the data and the best model ccdfs
[x y cdfvals] = CalcCDF_forKS(empdistr,char(models(1,:)),params{1,:},xmin,xmax);
x=reshape(x,numel(x),1); y=reshape(y,numel(y),1); cdfvals=reshape(cdfvals,numel(cdfvals),1);

ccdfs=zeros(length(x),N_plots+1);  %%initialize array
ccdfs(:,1:2)=[1-y,1-cdfvals];
figure('Name',datName); 
subplot(n1,n2,1)
hand(1) = loglog(x,ccdfs(:,1),'bo','MarkerSize',6,'MarkerFaceColor',[1 1 1]); hold on;
hand(2) = loglog(x,ccdfs(:,2),'k-','LineWidth',2); 
ylim([min(ccdfs(:,1)),max(ccdfs(:,1))]);
xlim([min(x),max(x)]);
title(char(models(1)),'Interpreter','None')
%%%if more than 1 model is supplied, plot them on same figure
if N_plots>1
    for j=2:N_plots
        subplot(n1,n2,j)
        hand(1) = loglog(x,ccdfs(:,1),'bo','MarkerSize',6,'MarkerFaceColor',[1 1 1]); hold on;
        [x y cdfvals] = CalcCDF_forKS(empdistr,char(models(j,:)),params{j,:},xmin,xmax);
        cdfvals=reshape(cdfvals,numel(cdfvals),1);
        ccdfs(:,j+1) =1-cdfvals;
        hand(j+1) = loglog(x,ccdfs(:,j+1),'k-','LineWidth',2); 
        ylim([min(ccdfs(:,1)),max(ccdfs(:,1))]);
        xlim([min(x),max(x)]);
        title(char(models(j)),'Interpreter','None')
        hold off
    end
end
% leg=legend(['Data'; models]);
% set(leg, 'Interpreter', 'none'); %%suppress LaTex interpretation of legend
% title(datName,'Interpreter', 'none');
% drawnow;
% hold off;
end

