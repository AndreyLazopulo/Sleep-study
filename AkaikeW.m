function [weights] = AkaikeW(AICs);
%The function takes in Akaike informatin criteria (AIC) for a number of 
%models and calculates the normalized Akaike weights for the models
%Sheyum 8.18.11
if iscell(AICs), AICs=cell2mat(AICs); end

deltas=AICs-min(AICs);
weights=exp(-deltas./2);%/(sum(exp(-deltas./2)));
weights=weights(:);
evratios=1./(weights./max(weights));
end

