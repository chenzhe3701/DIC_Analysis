% function [mask,hf] = plot_spline_mask(gb_dir, gb_s_pt, pt_pos, xMat, yMat)
% based on the inputs, plot grain boundary mask based on position matrices: xMat, yMat 
% chenzhe, 2018-05-15


function [mask,hf] = plot_spline_mask(gb_dir, gb_s_pt, pt_pos, xMat, yMat)

xy0 = xMat(1);   % determine the grid positions, e.g. [1 3 5 ...], or [11,16,21,...]
stepSize = yMat(2)-yMat(1);
mask = zeros(size(xMat));

hf = figure;
set(gca,'ydir','reverse');
hold on;
for ii=1:length(gb_s_pt)
    pos = pt_pos(gb_s_pt{ii},:);
    
    xmin = min(pos(:,1));
    xmax = max(pos(:,1));
    ymin = min(pos(:,2));
    ymax = max(pos(:,2));
%     xmin = xmin - mod(xmin-xy0,stepSize);
%     xmax = xmax - mod(xmax-xy0,stepSize);
%     ymin = ymin - mod(ymin-xy0,stepSize);
%     ymax = ymax - mod(ymax-xy0,stepSize);
    switch gb_dir{ii}
        case 'horizontal'
            pp = csapi(pos(:,1),pos(:,2));
            xp = xmin:stepSize/10:xmax;         % use smaller steps to ensure connectivity. chenzhe, 2018-06-14  
            yp = fnval(pp,xp);
        case 'vertical'
            pp = csapi(pos(:,2),pos(:,1));
            yp = ymin:stepSize/10:ymax;
            xp = fnval(pp,yp);
    end
    % make sure, round points of grain boundary to the grid
    xp = xp - mod(xp-xy0,stepSize);
    yp = yp - mod(yp-xy0,stepSize);
    plot(xp,yp,'--.k');
    
    % convert [x,y] to ind
    subx = (xp-xMat(1))/stepSize + 1;
    suby = (yp-yMat(1))/stepSize + 1;
    
    % force to be within boundary
    [nR,nC] = size(xMat);
    subx(subx<1) = 1;
    suby(suby<1) = 1;
    subx(subx>nC) = nC;
    suby(suby>nR) = nR;
    
    ind = sub2ind([nR,nC],suby,subx);
    mask(ind) = 1;
end

set(gca,'xlim',[xMat(1),xMat(end)],'ylim',[yMat(1),yMat(end)]);
end