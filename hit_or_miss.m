% I_logical = opening( I_input, K_kernel_default=ones(3) )
%
% chenzhe, 2018-05-18. Use standard name of this technique.


function I = hit_or_miss(I, C, D)

% make into 0 or 1
I(isnan(I)) = 0;

%I = double(erosion(I,C) & erosion(~I,D));
I(I~=0) = 1;
I = bwhitmiss(I, C, D);

end