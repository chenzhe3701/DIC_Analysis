% chenzhe 2018-02-28
% from 2 points, determine line equation y = kx + b
function [k, b] = pts_to_line_eq()
handleImline = imline(gca);

pos = wait(handleImline);

k = (pos(4)-pos(3))/(pos(2)-pos(1));
b = pos(3) - k*pos(1);

delete(handleImline);