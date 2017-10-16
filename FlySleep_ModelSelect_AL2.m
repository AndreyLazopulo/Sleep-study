function [funstatus] = FlySleep_ModelSelect_AL2(filename,data)
tic;
plotting=0;
binway=1;
funstatus=0;
% [f_name,p_name] = uigetfile({'*.TXT';'*.txt'},'Select files');
% filename=fullfile(p_name,f_name);
% data=dlmread(filename); % changed hdrload to dlmread
% indx=strfind(filename,'\');
% filename='data';
[~,f_name,~]=fileparts(filename);
disp(' ')
disp('Starting...')
disp(filename)
disp(' ')

inputarray=data(:,1);
inputarray(isnan(inputarray) | inputarray<=0)=[];  %remove NaN, 0 and neg values

distributdata = findandcut(inputarray,1);  %%%this function takes out events that occur only once
% anlbouts=inputarray;
minval=min(distributdata);
maxval=max(distributdata);
% fprintf('Original min bout: %g\n',minval)
% fprintf('Original max bout: %g\n',maxval)
% disp(' ')
disp(' ')
% if minval<0.5
%     anlbouts=anlbouts./minval;
%     fprintf('***!!!Scaling factor: %g***!!!\n',minval)
%     minval=min(anlbouts);
%     maxval=max(anlbouts);
%     fprintf('New min bout: %g\n',minval)
%     fprintf('New max bout: %g\n',maxval)  
%     disp(' ')
%     disp(' ')
% end

%plot data histogram
% figure;
% hold on
% h=histogram(log(anlbouts),0:0.2:10,'Normalization','probability');
% fitX=0.2/2:0.2:10;
% fitY=h.Values;
% plot(fitX,fitY,'o','Color','g')
% zerorem=inputdlg('how many zeros to remove?');
% zerorem=str2num(char(zerorem));

distributdata(distributdata<1)=[]; %%cut out out-of-range data

%make a cell array to store parameter values, sort them and do further
noofmodels=22; noofparameters=1;
n_rows=noofmodels+1; n_cols=8;
storecell = repmat({NaN},n_rows,n_cols);
storecell(1,:)={'Model','Params','LL'...
    'AICweight', 'LLtest p-value','LLtest R value','Gtest p-val','K-S p-val'};
mm=1; %model counter


%%3rd model: estimate parameters using the log-normal distr
mm=mm+1;
[mu,  sigma,  AIC_logno, ~, pval, LL]=Binned_LogNorm(distributdata,plotting,binway); % changed to my function
disp('Log-normal fit results: ')
fprintf('mu=%g ; sigma=%g \n',mu,sigma)
% fprintf('The p-value for log-normal fit is %g\n',pval)
fprintf('AIC=%g\n',AIC_logno)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'Log-Norm', [mu, sigma], LL, AIC_logno};
storecell(mm,7)={pval};
% 
% %6th model: truncated power-law
mm=mm+1;
[Alpha, AIC_tpl , ~, pval, LL] = Binned_TruncPL(distributdata,plotting,binway);
disp('Truncated Power-law fit results: ')
fprintf('alpha=%g\n',Alpha)
% fprintf('The p-value for truncated power-law fit is %g\n',pval)
fprintf('AIC=%g\n',AIC_tpl)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'PowerLaw_bound', Alpha, LL, AIC_tpl};
storecell(mm,7)={pval};
% 
% % %%%%the general Weibull distribution
mm=mm+1;
[alph, bet, lamd, AIC_gWei, pval, LL] = Binned_genWeibull(distributdata,plotting,binway); % changed to my function
disp('General Weibull fit results: ')
fprintf('alpha=%g ; beta=%g ; lamda=%g \n',alph,bet,lamd)
fprintf('AIC=%g\n',AIC_gWei)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'genWeibull', [alph, bet, lamd], LL, AIC_gWei};
storecell(mm,7)={pval};

%%7th model: power-law with exponential truncation
mm=mm+1;
[Aplambda_1, AIC_plexpo, ~, pval, LL] = Binned_PL1Expo(distributdata,plotting,binway);
disp('Power-law with Expo cutoff fit results: ')
fprintf('alpha=%g\n',Aplambda_1(1))
fprintf('lamda=%g\n',Aplambda_1(2))
%fprintf('The p-value for truncated power-law fit is %g\n',pval)
fprintf('AIC=%g\n',AIC_plexpo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'PLExpoCutff', Aplambda_1, LL, AIC_plexpo};
storecell(mm,7)={pval};

% %%2nd model: estimate parameters using the bounded exponential distr
mm=mm+1;
[mu, AIC_bex, ~, pval, LL] = Binned_1Expo_mle(distributdata,plotting,binway);
disp('Bounded Exponential fit results: ')
fprintf('mu=%g\n',mu)
% fprintf('The p-value for power-law fit is %g\n',pval)
fprintf('AIC=%g\n',AIC_bex)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'Expo_bound',mu,LL,AIC_bex};
storecell(mm,7)={pval};

%%***************
Nt=30; %%%******
% if fN==1 rngtwist = RandStream.create('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(rngtwist); end
%%%***************

%%2nd model: estimate parameters using the biexponential distrib
disp('Running model: 2 Exponentials')
disp(' ')
mm=mm+1;
[plambda_2, AIC_2Expo , ~, pval, LL] = Binned_2Expo_mle(distributdata,plotting,binway);
disp('Bi-Exponential fit results: ')
fprintf('w =%g,%g\n',plambda_2(1:2))
fprintf('lambda=%g,%g\n',plambda_2(3:4))
fprintf('AIC=%g\n',AIC_2Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'2Exponential', plambda_2, LL, AIC_2Expo};
storecell(mm,7)={pval};

%3rd model: estimate parameters using a 3 exponential distrib
disp('Running model: 3 Exponentials')
disp(' ')
mm=mm+1;
try
    [plambda_3, AIC_3Expo , ~, pval, LL] = Binned_3Expo_mle(distributdata,plotting,binway);
catch
    plambda_3=NaN(1,6);
    AIC_3Expo=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('3 Exponential fit results: ')
fprintf('w =%g,%g,%g\n',plambda_3(1:3))
fprintf('lambda=%g,%g,%g\n',plambda_3(4:6))
fprintf('AIC=%g\n',AIC_3Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'3Exponential', plambda_3, LL, AIC_3Expo};
storecell(mm,7)={pval};
% 
% %%4th model: estimate parameters using a 4 exponential distrib
disp('Running model: 4 Exponentials')
disp(' ')
mm=mm+1;
try
    [plambda_4, AIC_4Expo , ~, pval, LL] = Binned_4Expo_mle(distributdata,plotting,binway);
catch
    plambda_4=NaN(1,8);
    AIC_4Expo=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('4 Exponential fit results: ')
fprintf('w =%g,%g,%g,%g\n',plambda_4(1:4))
fprintf('lambda=%g,%g,%g,%g\n',plambda_4(5:8))
fprintf('AIC=%g\n',AIC_4Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'4Exponential', plambda_4, LL, AIC_4Expo};
storecell(mm,7)={pval};
% 
% %4th model: estimate parameters using a 5 exponential distrib
disp('Running model: 5 Exponentials')
disp(' ')
mm=mm+1;
try
    [plambda_5, AIC_5Expo , ~, pval, LL] = Binned_5Expo_mle(distributdata,plotting,binway);
catch
    plambda_5=NaN(1,10);
    AIC_5Expo=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('5 Exponential fit results: ')
fprintf('w =%g,%g,%g,%g,%g\n',plambda_5(1:5))
fprintf('lambda=%g,%g,%g,%g,%g\n',plambda_5(6:10))
fprintf('AIC=%g\n',AIC_5Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'5Exponential', plambda_5, LL, AIC_5Expo};
storecell(mm,7)={pval};

% %5th model: estimate parameters using a 6 exponential distrib mm=mm+1;
% [plambda_6, AIC_6Expo , ~, pval, LL] =Binned_6Expo_mle(distributdata,plotting,binway);
%  disp('6 Exponential fit results: ')
% fprintf('w =%g,%g,%g,%g,%g %g\n',plambda_6(1:6))
% fprintf('lambda=%g,%g,%g,%g,%g\n',plambda_6(7:12))
% fprintf('AIC=%g\n',AIC_6Expo) 
% fprintf('Gtest p-value=%g\n',pval)
% disp(' ')
% storecell(mm,1:4)={'6Exponential', plambda_6, LL, AIC_6Expo};
% storecell(mm,7)={pval};

%%%%*****STARTING THE POWER-LAW WITH SUM OF EXPONENTIAL SERIES*****

mm=mm+1;
disp('Running model:Power-law with 2 Exponentials')
disp(' ')
[Aplambda_2, AIC_Pl2Expo,~, pval, LL] = Binned_PL2Expo(distributdata,plotting,binway);
disp('Power-law with 2 Exponential results: ')
fprintf('alpha =%g\n',Aplambda_2(1))
fprintf('weights =%g,%g\n',Aplambda_2(2:3))
fprintf('lambda=%g,%g\n',Aplambda_2(4:5))
fprintf('AIC=%g\n',AIC_Pl2Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'PL2Exponential', Aplambda_2, LL, AIC_Pl2Expo};
storecell(mm,7)={pval};


mm=mm+1;
disp('Running model:Power-law with 3 Exponentials')
disp(' ')
try
    [Aplambda_3, AIC_Pl3Expo,~, pval, LL] = Binned_PL3Expo(distributdata,plotting,binway);
catch
    Aplambda_3=NaN(1,7);
    AIC_Pl3Expo=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Power-law with 3 Exponential results: ')
fprintf('alpha =%g\n',Aplambda_3(1))
fprintf('weights =%g,%g,%g\n',Aplambda_3(2:4))
fprintf('lambda=%g,%g,%g\n',Aplambda_3(5:7))
fprintf('AIC=%g\n',AIC_Pl3Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'PL3Exponential', Aplambda_3, LL, AIC_Pl3Expo};
storecell(mm,7)={pval};
% 
% 
mm=mm+1;
disp('Running model:Power-law with 4 Exponentials')
disp(' ')
try
    [Aplambda_4, AIC_Pl4Expo,~, pval, LL] = Binned_PL4Expo(distributdata,plotting,binway);
catch
    Aplambda_4=NaN(1,9);
    AIC_Pl4Expo=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Power-law with 4 Exponential results: ')
fprintf('alpha =%g\n',Aplambda_4(1))
fprintf('weights =%g,%g,%g,%g\n',Aplambda_4(2:5))
fprintf('lambda=%g,%g,%g,%g\n',Aplambda_4(6:9))
fprintf('AIC=%g\n',AIC_Pl4Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'PL4Exponential', Aplambda_4, LL, AIC_Pl4Expo};
storecell(mm,7)={pval};
% 
mm=mm+1;
disp('Running model:Power-law with 5 Exponentials')
disp(' ')
try
    [Aplambda_5, AIC_Pl5Expo, ~, pval, LL] = Binned_PL5Expo(distributdata,plotting,binway);
catch
    Aplambda_5=NaN(1,11);
    AIC_Pl5Expo=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Power-law with 5 Exponential results: ')
fprintf('alpha =%g\n',Aplambda_5(1))
fprintf('weights =%g,%g,%g,%g,%g\n',Aplambda_5(2:6))
fprintf('lambda=%g,%g,%g,%g,%g\n',Aplambda_5(7:11))
fprintf('AIC=%g\n',AIC_Pl5Expo)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'PL5Exponential', Aplambda_5, LL, AIC_Pl5Expo};
storecell(mm,7)={pval};

% mm=mm+1;
% disp('Running model:Power-law with 6 Exponentials')
% disp(' ')
% [Aplambda_6, AIC_Pl6Expo, ~, pval, LL] = Binned_PL6Expo(distributdata,plotting,binway);
% disp('Power-law with 6 Exponential results: ')
% fprintf('alpha =%g\n',Aplambda_6(1))
% fprintf('weights =%g,%g,%g,%g\n',Aplambda_6(2:7))
% fprintf('lambda=%g,%g,%g,%g\n',Aplambda_6(8:13))
% fprintf('AIC=%g\n',AIC_Pl6Expo)
% fprintf('Gtest p-value=%g\n',pval)
% disp(' ')
% storecell(mm,1:4)={'PL6Exponential', Aplambda_6, LL, AIC_Pl6Expo};
% storecell(mm,7)={pval};

% %%%%*****STARTING THE STREACHED EXPONENT WITH EXPONENTS SERIES *****
% 

mm=mm+1;
disp('Running model:Streached exponent')
disp(' ')
try
    [plambda_st, AIC_st , ~, pval, LL] = StrExpo_mle(distributdata,plotting,binway);
catch
    plambda_st=NaN(1,3);
    AIC_st=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Streached exponent results: ')
fprintf('alpha =%g\n',plambda_st(1))
fprintf('lambda=%g\n',plambda_st(3))
fprintf('AIC=%g\n',AIC_st)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'StrExpo', plambda_st(1:2:3), LL, AIC_st};
storecell(mm,7)={pval};


mm=mm+1;
disp('Running model:Streached exponent with 1 exponential')
disp(' ')
try
    [plambda_st1, AIC_st1 , ~, pval, LL] = StrExpo_1Expo(distributdata,plotting,binway);
catch
    plambda_st1=NaN(1,5);
    AIC_st1=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Streached exponent with 1 Exponential results: ')
fprintf('alpha =%g\n',plambda_st1(1))
fprintf('weights =%g,%g\n',plambda_st1(2:3))
fprintf('lambda=%g,%g\n',plambda_st1(4:5))
fprintf('AIC=%g\n',AIC_st1)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'Str1Exponential', plambda_st1, LL, AIC_st1};
storecell(mm,7)={pval};


mm=mm+1;
disp('Running model:Streached exponent with 2 exponentials')
disp(' ')
try
    [plambda_st2, AIC_st2 , ~, pval, LL] = StrExpo_2Expo(distributdata,plotting,binway);
catch
    plambda_st2=NaN(1,7);
    AIC_st2=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Streached exponent with 2 Exponential results: ')
fprintf('alpha =%g\n',plambda_st2(1))
fprintf('weights =%g,%g,%g\n',plambda_st2(2:4))
fprintf('lambda=%g,%g,%g\n',plambda_st2(5:7))
fprintf('AIC=%g\n',AIC_st2)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'Str2Exponential', plambda_st2, LL, AIC_st2};
storecell(mm,7)={pval};


mm=mm+1;
disp('Running model:Streached exponent with 3 exponentials')
disp(' ')
try
    [plambda_st3, AIC_st3 , ~, pval, LL] = StrExpo_3Expo(distributdata,plotting,binway);
catch
    plambda_st3=NaN(1,9);
    AIC_st3=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Streached exponent 3 Exponential results: ')
fprintf('alpha =%g\n',plambda_st3(1))
fprintf('weights =%g,%g,%g,%g\n',plambda_st3(2:5))
fprintf('lambda=%g,%g,%g,%g\n',plambda_st3(6:9))
fprintf('AIC=%g\n',AIC_st3)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'Str3Exponential', plambda_st3, LL, AIC_st3};
storecell(mm,7)={pval};

mm=mm+1;
disp('Running model:Streached exponent with 4 exponentials')
disp(' ')
try
    [plambda_st4, AIC_st4 , ~, pval, LL] = StrExpo_4Expo(distributdata,plotting,binway);
catch
    plambda_st4=NaN(1,11);
    AIC_st4=NaN;
    pval=NaN;
    LL=zeros(size(LL));
end
disp('Streached exponent with 4 Exponential results: ')
fprintf('alpha =%g\n',plambda_st4(1))
fprintf('weights =%g,%g,%g,%g,%g\n',plambda_st4(2:6))
fprintf('lambda=%g,%g,%g,%g,%g\n',plambda_st4(7:11))
fprintf('AIC=%g\n',AIC_st4)
fprintf('Gtest p-value=%g\n',pval)
disp(' ')
storecell(mm,1:4)={'Str4Exponential', plambda_st4, LL, AIC_st4};
storecell(mm,7)={pval};

% mm=mm+1;
% disp('Running model:Streached exponent with 5 exponentials')
% disp(' ')
% [plambda_st5, AIC_st5 , ~, pval, LL] = StrExpo_5Expo(distributdata,plotting,binway);
% disp('Streached exponent with 5 Exponential results: ')
% fprintf('alpha =%g\n',plambda_st5(1))
% fprintf('weights =%g,%g\n',plambda_st5(2:7))
% fprintf('lambda=%g,%g\n',plambda_st5(8:13))
% fprintf('AIC=%g\n',AIC_st5)
% fprintf('Gtest p-value=%g\n',pval)
% disp(' ')
% storecell(mm,1:4)={'Str5Exponential', plambda_st5, LL, AIC_st5};
% storecell(mm,7)={pval};
% 
% %%%%%%%*****************END OF MODEL TESTING************************
% %%%%%%%*****************END OF MODEL TESTING************************
% %%%%%%%*****************END OF MODEL TESTING************************
% 
%%Adjust size of storecell according to number of models tested
storecell(mm+1:end,:)=[];

%%sort the models according to their AICs, lowest at the top
storecell(2:end,:) = sortmodel(storecell(2:end,:),4,'ascend');

%%calculate the Akaike weights and replace AICs with the corresponding
%%weights
storecell(2:end,4) = num2cell(AkaikeW(storecell(2:end,4)));

%%do the log-likelihood test and calculate the p value. Using the model
%%with the highest Akaike weight as the 'best' one to test all the rest
storecell(find(isnan(cell2mat(storecell(2:end,4))))+1,:)=[];
nofloops=size(storecell,1)-1;
nest_mod=importdata('nested_models.csv');
mod_list=nest_mod.textdata(2:end,1);
for nn=1:nofloops
    
    [D KSp1] = KolmogSmirnov(distributdata,char(storecell(nn+1,1)),...
        storecell{nn+1,2},minval,maxval);
    
    KSD(nn,1)=D; KSp(nn,1)=KSp1;
end
for nnn=1:nofloops
    for nnnn=nnn:nofloops
        [R_val p_val] = LLratiotest(storecell{nnn+1,3},storecell{nnnn+1,3});
        model1=storecell(nnn+1,1); model2=storecell(nnnn+1,1);
        ind1=find(strcmp(mod_list,model1));
        ind2=find(strcmp(mod_list,model2));
        if nest_mod.data(ind1,ind2)==0
            R(nnn,nnnn)=R_val; p(nnn,nnnn)=p_val;
        elseif nest_mod.data(ind1,ind2)==1
            deg_fr=numel(storecell{nnn+1,2})-numel(storecell{nnnn+1,2});
            R(nnn,nnnn)=R_val; p(nnn,nnnn)=1-chi2cdf(abs(R_val),abs(deg_fr));
        end
    end
end
outputfile2=['LLtest_',f_name,'.csv'];
fID=fopen(outputfile2,'w');
fprintf(fID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\r\n','\',storecell{2:end,1});
for j=2:nofloops+1
fprintf(fID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\r\n',storecell{j,1},p(j-1,:));
end
fclose(fID);

storecell(2:end,5)=num2cell(p(:,1));
storecell(2:end,6)=num2cell(R(:,1));
% storecell(2:end,6)=num2cell(KSD(:,1));  %%update cell structure column
storecell(2:end,8)=num2cell(KSp(:,1));  %%update cell structure column

[Obs Expect] = GetObsExpect_AL(distributdata,char(storecell(2,1)),storecell{2,2},30);
n_params=numel(storecell{2,2});
str=char(storecell(2,1));
if strncmp(str(2:end), 'Exponential', 10),n_params=n_params-1;end;
if strncmp(str(4:end), 'Exponential', 10),n_params=n_params-1;end;
if strncmp(str(7:end), 'Exponential', 10),n_params=n_params-1;end;
[gof pval] = Gtest(Obs,Expect,n_params); %do the G-test

%%%calculate errors in parameter estimation via bootstrapping
%[paramerrors] = BootstrapFlyModels(anlbouts,char(storecell(2,1)),1000);

disp('Finishing...')
disp(filename)
fprintf('Max value of data is: %g\n',max(data(:,1)))
fprintf('Min value of data is: %g\n',min(data(:,1)))
fprintf('Censored values less than %g and more than %g\n',minval,maxval)
fprintf('Selected model: %s distribution, out of %g models\n',char(storecell(2,1)),mm-1)
fprintf('Parameter of the model:%g\n',storecell{2,2})
%fprintf('SE of parameter :%g\n',paramerrors)
fprintf('GofFit and p values are:%g and %g\n',gof,pval)

%%%%write to file all the text that appears on the screen
fid = fopen('report.txt', 'a');
fprintf(fid,'%s\n',' ');
fprintf(fid,'%s\n',filename);
fprintf(fid,'Number of data points: %f\n',numel(distributdata));
fprintf(fid,'Max value of data is: %f\n',max(data(:,1)));
fprintf(fid,'Min value of data is: %f\n',min(data(:,1)));
fprintf(fid,'Censored values less than %f and more than %f\n',minval,maxval);
fprintf(fid,'Selected model: %s distribution, out of %f models\n',char(storecell(2,1)),mm-1);
fprintf(fid,'Parameter of the model:%f\n',storecell{2,2});
%fprintf(fid,'SE of parameter :%g\n',paramerrors);
fprintf(fid,'GofFit and p values are:%f and %f\n',gof,pval);
fclose(fid);


%%%PLOTTING DATA AND FUNCTION%%%%%%
modelfunc = PlotDataModelCdf(distributdata, storecell(2:5,1),...
    storecell(2:5,2), minval,maxval,f_name,1);

%%%SAVE storecell TO FILE%%%
outputfile=['storecell_',f_name,'.csv'];
fID=fopen(outputfile,'w');
fprintf(fID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\r\n',storecell{1,1},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,2},storecell{1,3:4},storecell{1,5},storecell{1,5},storecell{1,5},storecell{1,5},storecell{1,7:end});
for j=2:nofloops
    fprintf(fID,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\r\n',storecell{j,1:2},NaN(13-length(storecell{j,2}),1),sum(storecell{j,3}),storecell{j,4},p(1:4,j-1)',storecell{j,7:end});
end
fclose(fID);

funstatus=1;
t=toc
end

