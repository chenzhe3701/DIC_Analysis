% hline = plot_spline_line(H,hv,stepSize)
% H{himpoint{1}, himpoint{2}, ...}
% hv = 'horizontal' or 'vertical'
% stepSize is step size
%
% chenzhe, 2018-06-13, modify to plot by 2 colors: black=horizontal,
% can use red=vertical, help to debug. Or change back to black

function hline = plot_spline_line(H,hv,stepSize)
    pos = H{1}.getPosition;
    switch hv
        case 'horizontal'
            xp = pos(1);
            for ii = 2:length(H)
                pos = [pos; H{ii}.getPosition];
                xp = [xp, xp(end)+stepSize : pos(ii,1)];
            end
            pp = csapi(pos(:,1),pos(:,2));
            yp = fnval(pp,xp);
            
            hline = plot(xp,yp,'--.k');
            
        case 'vertical'
            yp = pos(2);
            for ii = 2:length(H)
                pos = [pos; H{ii}.getPosition];
                yp = [yp, yp(end)+stepSize : pos(ii,2)];
            end
            pp = csapi(pos(:,2),pos(:,1));
            xp = fnval(pp,yp);
            
            hline = plot(xp,yp,'--.k');
    end
    
%     hline = plot(xp,yp,'--.k');
end