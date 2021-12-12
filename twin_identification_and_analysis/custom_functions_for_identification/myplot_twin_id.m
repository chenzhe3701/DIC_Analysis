
% stack a few plots together, and view interactively

function [] = myplot_twin_id(M,varargin)

p = inputParser;

addRequired(p,'M');
addOptional(p,'M2',zeros(size(M))*NaN);
addOptional(p,'M3',zeros(size(M))*NaN);
addOptional(p,'TF',zeros(size(M))*NaN);
addParameter(p,'X',repmat(1:size(M,2),size(M,1),1));
addParameter(p,'Y',repmat([1:size(M,1)]',1,size(M,2)));
addParameter(p,'alpha',1,@isnumeric);
addParameter(p,'handle',-1,@ishandle);
addParameter(p,'r',0);
parse(p,M,varargin{:});

M = p.Results.M;
M2 = p.Results.M2;
M3 = p.Results.M3;
TF = p.Results.TF;
X = p.Results.X;
Y = p.Results.Y;
f = p.Results.handle;
facealpha = p.Results.alpha;
ratio = p.Results.r;   % reduce ratio

if 0==ratio
    ratio = ceil(length(M)/3000);    % reduce ratio if necessary
end
if ratio > 1
    display(['matrix is big, use reduced ratio = ',num2str(ratio)]);
    M = M(1:ratio:end,1:ratio:end);
    M2 = M2(1:ratio:end,1:ratio:end);
    M3 = M3(1:ratio:end,1:ratio:end);
    TF = TF(1:ratio:end,1:ratio:end);
end

if ~isempty(M(isinf(M)))
    M(isinf(M)) = nan;
end
if ~isempty(M2(isinf(M2)))
    M2(isinf(M2)) = 0;
end
if ~isempty(M3(isinf(M3)))
    M3(isinf(M3)) = 0;
end


% M3(M3==0) = nan;
% M2(M2==0) = nan;

limit_x_low = min(X(:));
limit_x_high = max(X(:));
limit_y_low = min(Y(:));
limit_y_high = max(Y(:));

if f==-1
    f=figure;
end

clim = quantile(M(:),[0.005, 0.995]);
rmin = nanmin(M(:));

% Modify M to make grainboundary as -inf.
M(1==TF) = -inf;
M2(1==TF) = -inf;
M3(1==TF) = -inf;

% modify the colormap to make smallest black, biggest white.
colormap([0 0 0; colormap; 1 1 1]);


hold on;
set(f,'position',[50,50,800,600]);
% title(strrep(inputname(1),'_','\_'));

ax1 = axes;
% Let the alpha of nan points to be 0.
Malpha = ones(size(M));
Malpha(isnan(M)) = 0;
imagesc([X(1),X(end)],[Y(1),Y(end)],M,'alphadata',Malpha);
try
    caxis(gca, [clim(1), clim(2)]);
end


ax2 = axes;
M2alpha = ones(size(M2));
M2alpha(isnan(M2)) = 0;
imagesc([X(1),X(end)],[Y(1),Y(end)],M2,'alphadata',M2alpha);
try
    caxis(gca, [-0.1, 6.1]);
end


ax3 = axes;
M3alpha = ones(size(M3));
M3alpha(isnan(M3)) = 0;
imagesc([X(1),X(end)],[Y(1),Y(end)],M3,'alphadata',M3alpha);
try
    caxis(gca, [0, 5]);
end


% the following handles slider. Can turn off
s = uicontrol('Parent',f,...
    'Style','slider',...
    'Units','normalized',...
    'Position',[0.93, 0.15, 0.03, 0.45],...
    'Min',1,'Max',3,'Value',3,'SliderStep',[1/2, 1/2]);

s.Callback = {@callback_slider,ax1,ax2,ax3};

b = uicontrol('Parent',f,'Style','pushbutton','String','select grain','Units','normalized','Position',[0.93, 0.7, 0.06, 0.05]);
b.Callback = {@callback_pushbutton_1};

b2 = uicontrol('Parent',f,'Style','pushbutton','String','input correction','Units','normalized','Position',[0.93, 0.78, 0.06, 0.05]);
b2.Callback = {@callback_pushbutton_2};

end

function callback_slider(source,event,ax1,ax2,ax3)
% get value of the slider
value = round(source.Value);

switch value
    case 1
        axes(ax1);
    case 2
        axes(ax2);
    case 3
        axes(ax3);
end    
    
end

function callback_pushbutton_1(source,event)
    evalin('base','run twin_id_script_select_grain.m;');
end

function callback_pushbutton_2(source,event)
    evalin('base','run twin_id_script_input_correction.m;');
end