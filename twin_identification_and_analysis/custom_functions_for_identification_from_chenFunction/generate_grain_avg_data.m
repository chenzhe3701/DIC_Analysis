% [avgOut, stdOut] = generate_grain_avg_data(ID,gID,inputData,goodValueThreshold,sigma)
%
% ID, inputData are matrices
% gID is a vector containing IDs of the interested grains
% if grain has >goodValueThreshold (a pct) of points with valid value, then calculate average values for this grain
% Otherwise, the value for this grain is NaN.
%
% Zhe Chen, 2015-08-12 revised.

function [avgOut, stdOut] = generate_grain_avg_data(ID,gID,inputData,goodValueThreshold,sigma)
if nargin < 5
    sigma = zeros(size(ID));
    disp('no sigma value provided, assume all data points valid');
end
avgOut = zeros(length(gID),1)*nan;
stdOut = zeros(length(gID),1)*nan;
for iGrain = 1:length(gID)
    currentID = gID(iGrain);
    indexTemp = (ID == currentID);
    indexTemp2 = indexTemp&(sigma~=-1);
    if sum(indexTemp2(:))/sum(indexTemp(:)) > goodValueThreshold   % if most of the grain has valid strain measurement
        valueTemp = inputData(indexTemp2);
        avgOut(iGrain) = nanmean(valueTemp(:));
        stdOut(iGrain) = nanstd(valueTemp(:));
    end
end

end