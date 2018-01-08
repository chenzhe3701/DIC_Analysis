% use sgolayfilt twice to perform filtering in 2d
%
% chenzhe, 2018-01-06

function y=sgolayfilt2(x,order,framelen)

y = sgolayfilt(sgolayfilt(x,order,framelen,[],1),order,framelen,[],2);