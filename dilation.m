% I = dilation( I_input, K_kernel_default=ones(3) )
%
% Use imdilate, but take care of NaNs (and -inf ?)
% chenzhe, 2018-05-18. Use standard name of this technique.


function I = dilation(I, K)

% can specify a kernel
if ~exist('K','var')
    K = ones(3);
end

% make into 0 or 1
I(isnan(I)) = 0;

% Is this needed?
% I(I~=0) = 1;

I = imdilate(I,K);

% Is this needed? Maybe not.
% I(isinf(I)) = 1;

end