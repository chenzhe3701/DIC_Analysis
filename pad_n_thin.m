% [I] = thin_once(I,nTimes) 
% For special need. First padd I with extra boundary of '0's, then thin,
% return the same size .
% Actually it can thin 'nTimes'.
% chenzhe, 2018-05-20.
function [I] = pad_n_thin(I,nTimes)

if ~exist('nTimes','var')
    nTimes = 1;
end

% make into 0 or 1
I(isnan(I)) = 0;
I(I~=0) = 1;

I = padarray(I,[1,1],'pre');
I = padarray(I,[1,1],'post');
I = bwmorph(I,'thin',nTimes);
I = I(2:end-1,2:end-1);
end