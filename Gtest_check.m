function all_pval=Gtest_check
anlbouts=dlmread('w1118_all.txt');
anlbouts = findandcut(anlbouts,1);
ff=importdata('storecell_w1118_all.csv');
models=ff.textdata(2:end,1);
allparams=ff.data(1:end,1:13);
all_pval=zeros(length(models),10);
for k=1:20
    for i=1:length(models)
        params=allparams(i,~isnan(allparams(i,:)));
        [Obs Expect] = GetObsExpect_AL2(anlbouts,char(models(i)),params,k);
        n_params=numel(params);
        str=char(models(i));
        if strncmp(str(2:end), 'Exponential', 10),n_params=n_params-1;end;
        if strncmp(str(4:end), 'Exponential', 10),n_params=n_params-1;end;
        if strncmp(str(7:end), 'Exponential', 10),n_params=n_params-1;end;
        [gof pval] = Gtest(Obs,Expect,n_params); %do the G-test
        if ~isempty(pval)
            all_pval(i,k)=pval;
        end
    end
end
for j=1:length(models)
    figure;
    plot(all_pval(j,:));
    title(char(models(j)),'Interpreter','none')
    xlabel('smallest bin (min)')
    ylabel('p-value')
end
    