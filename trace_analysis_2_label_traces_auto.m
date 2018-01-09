% Start with analyzing more data, such as quantiles, ...

clear;
addChenFunction;
dicPath = uigetdir('','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);

gIDwithTrace = gID(~isnan(gExx));
% gIDwithTrace = [89,129,135,191,201,210,327,401,422,553];        % WE43 T5 #7, new
% gIDwithTrace = [302,151,186,191,296,431,572,1211];            % WE43_T6_C1

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 1;   % elongation levels to analyze. 0-based.
iE_stop = 7;
resReduceRatio = 3;         % to save space, reduce map resolution
grow_boundary_TF = 0;       % whether to grow boundary to make it thicker
% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'

notes = struct('atEdge',[],'likeTwin',[],'tooSmall',[]);
% This is a special note for data processing, for specific samples, to give grain labels
% useful fields will be 'tooSmall', 'atEdge', 'likeTwin'
% notes = load('T5#7_traceAnalysisNotes');
% 2016-08-15. create a field called 'noTrace'. noTrace = 1 = too small. noTrace = 2 = at edge .

% end of modify settings part 1 ------------------------------------------------------------------------------------------------------------------------------------
save([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat'],...
    'dicPath','dicFiles',...
    'STOP','iE_start','iE_stop','resReduceRatio','grow_boundary_TF','f1','f2','neighbor_elim','twinTF_text','notes','gIDwithTrace',...
    '-append');

%%
iE = 5;
strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);            % change the prefix and name of DIC file ----------------------------------------------------

load(strainFile,'exx','exy','eyy');     % Look at exx, but this can be changed in the future.   % ----------------------------------------------------------------------------------

%%
for iS = 1%:length(gIDwithTrace)
    close all;
    ID_current = gIDwithTrace(iS);              % id of current grain
    ID_current = 89;
    
    ind_current = find(ID_current == gID);    % an index of row
    phi1_current = gPhi1(ind_current);
    phi_current = gPhi(ind_current);
    phi2_current = gPhi2(ind_current);
    
    ID_neighbor = gNeighbors(ind_current,:);
    ID_neighbor = ID_neighbor((ID_neighbor~=0)&(ID_neighbor~=neighbor_elim));
    
    % find index range of a small matrix containing the grain of interest
    ind_pool = ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_pool, 1), 1, 'first');
    indC_max = find(sum(ind_pool, 1), 1, 'last');
    indR_min = find(sum(ind_pool, 2), 1, 'first');
    indR_max = find(sum(ind_pool, 2), 1, 'last');
    
    nRow = indR_max - indR_min + 1;
    nColumn = indC_max - indC_min + 1;
    
    e_current = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later --------------------
    boundaryTF_current = boundaryTF(indR_min:indR_max, indC_min:indC_max);
    x_current = X(indR_min:indR_max, indC_min:indC_max);
    y_current = Y(indR_min:indR_max, indC_min:indC_max);
    ID_map_current = ID(indR_min:indR_max, indC_min:indC_max);
    
    e_grain = e_current;
    e_grain(ID_map_current~=ID_current) = 0;  % 'e_grain' is strain of This grain. 'e_current' is strian of this region.
    e_grain(isnan(e_grain)) = 0;
    
    % calculate Schmid factor.
    % sf_mat = [#, SF, angle_XtoY, trace_x_end, trace_y_end]
    % burgersXY = [burgers_X, burgers_Y, ratio].
    % ---------------------------------------- select proper setting for analysis, such as material, twin, stress ---------------------------------------------------------
    
    [sf_mat, sf_mat_sorted, burgersXY] = trace_analysis_TiMgAl([phi1_current,phi_current,phi2_current],[-90,180,0],[0,0,0],stressTensor,sampleMaterial,twinTF_text);
    
    
   
    figureHandle_1 = myplot(x_current,y_current,e_grain);
    
    % hough transform, houghpeaks the peaks, houghlines the line segments
    [h,t,r] = hough(e_grain);
    [tg,rg] = meshgrid(t,r);
    % figureHandle_2 = myplot(tg,rg,h);
    figureHandle_2 = figure;
    imshow(h,[],'XData',t,'YData',r,'colormap',parula); axis on; axis square;
    
    peaks = houghpeaks(h,10);   %,'NHoodSize',[floor(size(h,1)/2)*2+1,3]);
    set(0,'currentfigure',figureHandle_2); hold on;
    for k = 1:size(peaks,1)
        xy = [tg(peaks(k,1),peaks(k,2)),rg(peaks(k,1),peaks(k,2))];
        plot3(xy(1),xy(2),max(h(:)),'s','LineWidth',1,'Color','k')
    end
    
    lines = houghlines(e_grain,t,r,peaks);
    % lines(k).point1/2 is in fact the index (index_c, index_r)
    % show the extracted lines
    set(0,'currentfigure',figureHandle_1); hold on;
    for k = 1:length(lines)
        xy = [x_current(lines(k).point1(2),lines(k).point1(1)),y_current(lines(k).point1(2),lines(k).point1(1));...
            x_current(lines(k).point2(2),lines(k).point2(1)),y_current(lines(k).point2(2),lines(k).point2(1))];
        plot3(xy(:,1),xy(:,2),1*ones(size(xy,1),1),'LineWidth',2,'Color','green')
    end
    
    
end






















