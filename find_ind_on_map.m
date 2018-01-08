% with a myplot() plotted map of strain fields/grain boundary overlays, ...
% manually crop an area, and return the [indr[], indc[]] which can be used
% later for croping data
%
% chenzhe, 2018-01-06

function [indR,indC] = find_ind_on_map(X,Y)
h = imrect;
position = wait(h)
delete(h);
[~,indr1] = min(abs(Y(:,1)-position(2)));
[~,indr2] = min(abs(Y(:,1)-position(2)-position(4)));
[~,indc1] = min(abs(X(1,:)-position(1)));
[~,indc2] = min(abs(X(1,:)-position(1)-position(3)));
indR = indr1:indr2;
indC = indc1:indc2;
end