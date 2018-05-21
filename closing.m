% I = closing( I_input, K_kernel_default=ones(3) )
%
% chenzhe, 2018-05-18. Use standard name of this technique.


function I = closing(I, K)

% can specify a kernel
if ~exist('K','var')
    K = ones(3);
end

% make into 0 or 1
I(isnan(I)) = 0;
% I(I~=0) = 1;

% I = erosion(dilation(I,K), K);
I = imclose(I, K);

end