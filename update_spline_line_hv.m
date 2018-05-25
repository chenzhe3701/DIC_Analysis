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
        
        % In case the control points move from larger to smaller value ...
        p2 = himpoint{end}.getPosition;
        x_p2 = p2(1);
        if x_p2 < xp
            stepSize = -stepSize;
        end
        
        for ii = 2:length(himpoint)
            pos = [pos; himpoint{ii}.getPosition];
            xp = [xp, xp(end)+stepSize:sign(stepSize):pos(ii,1)];
        end
        pp = csapi(pos(:,1),pos(:,2));
        yp = fnval(pp,xp);
    case 'vertical'
        yp = pos(2);
        
        % In case the control points move from larger to smaller value ...
        p2 = himpoint{end}.getPosition;
        y_p2 = p2(2);
        if y_p2 < yp
            stepSize = -stepSize;
        end
        
        for ii = 2:length(himpoint)
           pos = [pos; himpoint{ii}.getPosition];
           yp = [yp, yp(end)+stepSize:sign(stepSize):pos(ii,2)];
        end
        pp = csapi(pos(:,2),pos(:,1));
        xp = fnval(pp,yp);
end


hline.XData = xp;
hline.YData = yp;

% update if necessary. But may be slow. Can comment this.
pt_pos = evalin('base','pt_pos;');
h = evalin('base','h;');
for ii = 1:size(pt_pos,1)
    try
        pt_pos(ii,:) = h{ii}.getPosition;
    end
end
assignin('base','pt_pos',pt_pos);
assignin('base','h',h);

end