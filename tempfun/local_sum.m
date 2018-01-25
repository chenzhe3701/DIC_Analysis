%
% imagine we have a 2d-cumsum [AB;CD], then D = ABCD-AC-AB+A
% [AB;CD] is 'fen kuai ju zhen'
% 
% comparing to conv2, this is much faster.  But does not work well with nan
% chenzhe, 2018-01-23, based on normxcorr2A

function s = local_sum(A,m,n)
B = A;
B = padarray(B,[m,n],'pre');
B = padarray(B,[m-1,n-1],'post');
B = cumsum(B,1);
B = cumsum(B,2);
s = B(1+m:end,1+n:end) - B(1+m:end, 1:end-n) - B(1:end-m, 1+n:end) + B(1:end-m, 1:end-n);
end