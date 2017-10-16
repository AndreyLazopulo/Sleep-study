function [trimmed] = findandcut(inputarray,freq )
%%%This function takes out all elements that occur as often as or less often
%%%than 'freq' in an 'inputarray'. Default value of 'freq' is 1.

if isempty(freq), freq=1; end;

trimmed=inputarray;  %initialize
[cutvalues numvalues]=count_unique(inputarray);

%construct vector that contain elements
%that occur only 'freq' or fewer no. of times
cutvalues(numvalues>freq)=[];  

%%%%cut those elements out
if numel(cutvalues)>0
    for i=1:numel(cutvalues)
        trimmed(trimmed==cutvalues(i))=[];
    end
end

end

