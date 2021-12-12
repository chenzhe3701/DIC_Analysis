%% [***] script to change to double axis.
changeColor = 0;

if changeColor==0
    %  (1) move labels of axis-1 to top and right
    a = gca; title(a,'');
    set(a, 'fontsize',24);
    set(a,'XAxisLocation','top','YAxisLocation','right', 'Color','none');
    xlabel(a,'X (mm)','color','k');
    ylabel(a,'Y (mm)','color','k');
    
    % record original tick
    xTick = a.XTick;
    yTick = a.YTick;
    xTickLabel = a.XTickLabel;
    yTickLabel = a.YTickLabel;
    
    % artificially change tick label of axis-1
    scale = 4096/0.36;
    set(a,'XTick',[0:5]*scale,'XTickLabel',{'0','1','2','3','4','5'},'YTick',[0:3]*scale,'YTickLabel',{'0','1','2','3'});
    
    % (2) make axis-2, and make the labels back to the left and bottom
    pause(.1);
    ax2 = copyobj(a,a.Parent);
    ax2.ActivePositionProperty = 'Position';
    set(ax2, 'XAxisLocation','bottom','YAxisLocation','left', 'Color','none');
    xlabel(ax2,'X (\times10^4 pixels)','color','k');
    ylabel(ax2,'Y (\times10^4 pixels)','color','k');
    
    set(ax2,'XTick',xTick,'YTick',yTick);
    set(ax2,'XTickLabel', xTickLabel);
    set(ax2,'YTickLabel', yTickLabel);
    title(ax2,'');
else
    %  (1) move labels of axis-1 to top and right
    a = gca; title(a,'');
    set(a, 'fontsize',24);
    set(a,'XAxisLocation','top','YAxisLocation','right', 'Color','none');
    xlabel(a,'X (mm)','color','b');
    ylabel(a,'Y (mm)','color','b');
    
    % record original tick
    xTick = a.XTick;
    yTick = a.YTick;
    xTickLabel = a.XTickLabel;
    yTickLabel = a.YTickLabel;
    
    % artificially change tick label of axis-1
    scale = 4096/0.36;
    set(a,'XTick',[0:5]*scale,'XTickLabel',{'0','1','2','3','4','5'},'YTick',[0:3]*scale,'YTickLabel',{'0','1','2','3'});
    
    labels = [];
    for ii = 1:length(a.XTickLabel)
        labels{ii} = ['\color{blue}',num2str(a.XTickLabel{ii})];
    end
    set(a,'XTickLabel',labels);
    
    labels = [];
    for ii = 1:length(a.YTickLabel)
        labels{ii} = ['\color{blue}',num2str(a.YTickLabel{ii})];
    end
    set(a,'YTickLabel',labels);
    
    % (2) make axis-2, and make the labels back to the left and bottom
    pause(.1);
    ax2 = copyobj(a,a.Parent);
    ax2.ActivePositionProperty = 'Position';
    set(ax2, 'XAxisLocation','bottom','YAxisLocation','left', 'Color','none');
    xlabel(ax2,'X (\times10^4 pixels)','color','k');
    ylabel(ax2,'Y (\times10^4 pixels)','color','k');
    
    set(ax2,'XTick',xTick,'YTick',yTick);
    labels = [];
    for ii = 1:length(xTickLabel)
        labels{ii} = ['\color{black}',num2str(xTickLabel{ii})];
    end
    set(ax2,'XTickLabel',labels);
    
    labels = [];
    for ii = 1:length(yTickLabel)
        labels{ii} = ['\color{black}',num2str(yTickLabel{ii})];
    end
    set(ax2,'YTickLabel',labels);
    
    title(ax2,'');
end
