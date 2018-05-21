% I = erosion( I_input, K_rotated_kernel_default=ones(3) )
%
% Use imerode, but take care of NaNs (and inf ?)
% Note that imerode use
% chenzhe, 2018-05-18. Use standard name of this technique.


function I = erosion(I, K)

% can specify a kernel
if ~exist('K','var')
    K = ones(3);
end

% make into 0 or 1
I(isnan(I)) = 0;

% Is this needed?
% I(I~=0) = 1;

I = imerode(I,K);

% Is this needed? Maybe not.
% I(isinf(I)) = 1;

end