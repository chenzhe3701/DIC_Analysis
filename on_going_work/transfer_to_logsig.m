% map x into X, i.e., X=f(x), so that logsig(f(x)) is s-shaped probability
% logsig(f(xm)) = 0.5
% logsig(f(xp)) = logsig(xd) = p
% chenzhe, 2018-09-23
%
% Explanation:
%
% logsig = 1 / 1+exp(-x) = p
% 
% to make [xp] corresponding to a probability of [p]
% [xm] corresponding to a probability of  0.5
% 
% (0) We have to solve for xd that correspond to a probability of p
% 1 / 1+exp(-xd) = p --> xd = ln(p/1-p)
% --> this xd is called 'logit' or 'log-odds'
%
% (1) we need first shift to re-zero:
% X1 = x - xm
% 
% (2) then scale
% X2 = F(x-xm), where F = (xd-0)/(xp-xm)
% 
% So, the whole transfer should be
% 
% X = (xd-0) / (xp-xm) * (x-xm)
%   = ln(p/(1-p)) / (xp-xm) * (x-xm)


function [x,xd] = transfer_to_logsig(x, xm, xp, p)

xd = nan;

if xp==xm
    warning('xp should not be equal to xm');
    return;
end

xd = log(p/(1-p));
x = xd / (xp-xm) * (x-xm);
