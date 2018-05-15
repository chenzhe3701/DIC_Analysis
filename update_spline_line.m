% The input 'p' is the output from

function [] = update_spline_line(hline, himpoint, stepSize)

pos = himpoint{1}.getPosition;
xp = pos(1);
for ii = 2:length(himpoint)
    pos = [pos; himpoint{ii}.getPosition];
    xp = [xp, linspace(xp(end)+stepSize, pos(ii,1), 100)];
end
pp = csapi(pos(:,1),pos(:,2));


yp = fnval(pp,xp);

hline.XData = xp;
hline.YData = yp;

end