% The input should include 'p', which is the output an 'himpoint' when
% position is update. However, all related point handles are input in
% 'himpoint', so 'p' is ignored.
%
% hline: handle of a line by 'plot'
% himpoint: cell array of impoint handles.
% v: whether 'horizontal' fit or 'vertical' fit
% stepSize: distance between points

function [] = update_spline_line_hv(hline, himpoint, hv, stepSize)

pos = himpoint{1}.getPosition;

switch hv
    case 'horizontal'
        xp = pos(1);
        for ii = 2:length(himpoint)
            pos = [pos; himpoint{ii}.getPosition];
            xp = [xp, xp(end)+stepSize : pos(ii,1)];
        end
        pp = csapi(pos(:,1),pos(:,2));
        yp = fnval(pp,xp);
    case 'vertical'
        yp = pos(2);
        for ii = 2:length(himpoint)
           pos = [pos; himpoint{ii}.getPosition];
           yp = [yp, yp(end)+stepSize : pos(ii,2)];
        end
        pp = csapi(pos(:,2),pos(:,1));
        xp = fnval(pp,yp);
end


hline.XData = xp;
hline.YData = yp;

end