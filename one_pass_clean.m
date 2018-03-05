% chenzhe, 2018-03-04
% use one_pass_label to generate unique ID map of input matrix
% for the uniqueID map, disable the ones with grain size smaller than the threshold
% then return the valid part (i.e., with grain size > threshold) of the input matrix
% min_size, default = 0.25% of input data size

function cleaned = one_pass_clean(raw, min_size)

if ~exist('min_size','var')
    min_size = sum(raw(:)>0)*0.0025;
end


label = one_pass_label(raw);
uniqueLabel = unique(label(:));
uniqueLabel(isnan(uniqueLabel)) = [];
gSize = [];
for ii = 1:length(uniqueLabel)
    gSize(ii) = sum(label(:)==uniqueLabel(ii));
end
validLabel = uniqueLabel(gSize>min_size);
validMap = zeros(size(raw));
validMap(ismember(label,validLabel(:))) = 1;
cleaned = raw .* validMap;
end