function average_sleep=aversleep(interval,winsize,wrapinter)

[filename,pathname]=uigetfile('.txt', 'Select data files','MultiSelect','on');
    % find number of selected files
    if iscell(filename)
        NumberOfFiles=length(filename);
    elseif filename ~= 0
        NumberOfFiles = 1;
    else
        NumberOfFiles = 0;
    end
    average_sleep=zeros(5000000,1);
for i=1:NumberOfFiles
    
    if NumberOfFiles == 1
        file1=fullfile(pathname,filename);
    else
        file1=fullfile(pathname,filename(i));
        file1=char(file1);
    end
    
    % open data file
    [~,X]=hdrload(file1);
    X=X(1:end);
    if isempty(X)==1
        continue
    end
    S=sleeptrace(X,5/interval,0);
    if wrapinter=='nowrap'
        new_sleep=RunningSleep2(S,interval,winsize,length(S)/(60*interval));
        N=min(length(average_sleep),length(new_sleep));
        average_sleep=average_sleep(1:N)+new_sleep(1:N);
    else
        average_sleep=average_sleep+RunningSleep2(S,interval,winsize,wrapinter);
    end
end
average_sleep=average_sleep./NumberOfFiles;