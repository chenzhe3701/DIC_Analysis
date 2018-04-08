% script
% chenzhe, 20180308, prepare figures for TMS

addChenFunction;
grainDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\Grain_1144_data_for_paper_ppt','Folder that contains the grain data'),'\'];

dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

% load previous data and settings
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis'],'X','Y','boundaryTF','boundaryTFB','ID','gID','gExx','exx');
% load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
gIDwithTrace = gID(~isnan(gExx));


% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

% file name prefixes
f1 = 'WE43_T6_C1_s';
f2 = '_';

[ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
addChenFunction;

%% Necessary information for cNumMaps etc.
egc = 211441;   % [iE,grainID,iCluster] format
iC_target = mod(egc,10);
egc = (egc-iC_target)/10;
ID_target = mod(egc,10000);
iE_target = (egc-ID_target)/10000;
ID_current=ID_target;

for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'clusterNumMap','stru','clusterNumMapCleaned');
    cNumMaps{iE} = clusterNumMap;
    cNumMaps_cleaned{iE} = clusterNumMapCleaned;
    struCell{iE} = stru;
end

iS = find(arrayfun(@(x) x.gID == ID_target,struCell{iE_start}));    % just find iS from any of the stru
% a good method to, from the starting iE and iC, fill all useful iEs and iCs
iE_list = iE_target;
iC_list = iC_target;
while 0 ~= struCell{iE_list(1)}(iS).preCluster(iC_list(1))
    iC_list = [struCell{iE_list(1)}(iS).preCluster(iC_list(1)), iC_list];
    iE_list = [iE_list(1)-1, iE_list];
end

while 0 ~= struCell{iE_list(end)}(iS).postCluster(iC_list(end))
    iC_list = [iC_list,struCell{iE_list(end)}(iS).postCluster(iC_list(end))];
    iE_list = [iE_list, iE_list(end)+1];
end





%% [1] plot the strain, RGB map at the strain level you are interested in: iE = 5

% iE of interest ----------------------------- !!! -------------------------------------
iE = 5;

fname = [f1,'_all_grain_',num2str(ID_current),'_local_map.mat'];
load(fullfile(grainDataPath,fname),'data');
s = data(iE);

indR_min = s.indR_min;
indR_max = s.indR_max;
indC_min = s.indC_min;
indC_max = s.indC_max;
boundaryTF_local = s.boundaryTF_local;
x_local = s.x_local;
y_local = s.y_local;
ID_current = s.ID_current;
ID_local = s.ID_local;
exx_local = s.exx_local;
exy_local = s.exy_local;
eyy_local = s.eyy_local;
sigma_local = s.sigma_local;
clusterNumMapLocal = s.clusterNumMapLocal;
clusterNumMapCleanedLocal = s.clusterNumMapCleanedLocal;

%% (1) plot/adjust strain map here -------------------------------------------------
[f,a,c] = myplot(x_local-x_local(1), y_local-y_local(1), exx_local); % caxis([-0.1, 0.00]);% caxis([-0.14, 0.07]);
title(a,'\epsilon_x_x');
set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
axis equal;
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_exx.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);

[f,a,c]=myplot(x_local-x_local(1), y_local-y_local(1), exy_local); % caxis([-0.03, 0.03]);% caxis([-0.07, 0.07]);
title(a,'\epsilon_x_y');
set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
axis equal;
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_exy.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);

[f,a,c]=myplot(x_local-x_local(1), y_local-y_local(1), eyy_local); % caxis([-0.02, 0.06])% caxis([-0.07, 0.14]);
title(a,'\epsilon_y_y');
set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
axis equal;
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eyy.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);













%% (2) convert to RGB image
img_size = 227;
[nR,nC] = size(clusterNumMapLocal);
sz = max(nR,nC);
img_local_cNum = zeros(sz)*nan;
img_exx_local = zeros(sz)*nan;
img_exy_local = zeros(sz)*nan;
img_eyy_local = zeros(sz)*nan;
if nR<nC
    img_local_cNum(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = clusterNumMapLocal;
    img_exx_local(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = exx_local;
    img_exy_local(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = exy_local;
    img_eyy_local(round((sz-nR)/2)+1:round((sz-nR)/2)+nR, :) = eyy_local;
else
    img_local_cNum(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = clusterNumMapLocal;
    img_exx_local(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = exx_local;
    img_exy_local(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = exy_local;
    img_eyy_local(:, round((sz-nC)/2)+1:round((sz-nC)/2)+nC) = eyy_local;
end

ratio = sz/img_size;
img_local_cNum = img_local_cNum(ratio:ratio:end,ratio:ratio:end);
img_exx_local = img_exx_local(ratio:ratio:end,ratio:ratio:end);
img_exy_local = img_exy_local(ratio:ratio:end,ratio:ratio:end);
img_eyy_local = img_eyy_local(ratio:ratio:end,ratio:ratio:end);

% scale the value, prepare for image. !!!
img_exx_local = mat_to_image(img_exx_local, [-0.14, 0.07], 'index');
img_exy_local = mat_to_image(img_exy_local, [-0.07, 0.07], 'index');
img_eyy_local = mat_to_image(img_eyy_local, [-0.07, 0.14], 'index');

%% do the RGB image of here ---------------------------------------------------------
% [A] cluster by cluster
stru = struCell{iE};
nCluster = length(stru(iS).cLabel);
for iCluster = 1:nCluster
    cNum = stru(iS).cLabel(iCluster);
    % imgLocal = uint8((img_local_cNum==cNum)*255);
    exx_R_Local = uint8(  img_exx_local.*(img_local_cNum==cNum) *255);
    exy_G_Local = uint8(  img_exy_local.*(img_local_cNum==cNum) *255);
    eyy_B_Local = uint8(  img_eyy_local.*(img_local_cNum==cNum) *255);
    
    bkgd = uint8(zeros(size(exx_R_Local)));
    imgR = cat(3,exx_R_Local, bkgd, bkgd);
    imgG = cat(3,bkgd, exy_G_Local, bkgd);
    imgB = cat(3,bkgd, bkgd, eyy_B_Local);
    % imgLocal = cat(3,imgLocal,imgLocal,imgLocal);
    imgLocal = cat(3,exx_R_Local,exy_G_Local,eyy_B_Local);
    
    % (1) R channel
    f = figure;
    ax = axes;
    colorMap = zeros(3,64);
    colorMap(1,:) = linspace(0,1,64);
    colorMap = colorMap';
    image(ax,imgR);
    title(ax,'R channel');
    colormap(colorMap);
    c=colorbar;
    title(c,'\epsilon_x_x')
    c.Ticks = [0,0.33,0.66,1];
    set(c,'TickLabels',{'-0.14','-0.07','0','0.07'});
    set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
    axis equal; axis off;
    
    imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_c',num2str(iCluster),'_R','.tif']);
    print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    close(f);
    
    % (2) G channel
    f = figure;
    ax = axes;
    colorMap = zeros(3,64);
    colorMap(2,:) = linspace(0,1,64);
    colorMap = colorMap';
    image(ax,imgG);
    title(ax,'G channel');
    colormap(colorMap);
    c=colorbar;
    title(c,'\epsilon_x_y');
    c.Ticks = [0,0.5,1];
    set(c,'TickLabels',{'-0.07','0','0.07'});
    set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
    axis equal; axis off;
    
    imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_c',num2str(iCluster),'_G','.tif']);
    print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    close(f);
    
    % (3) B channel
    f = figure;
    ax = axes;
    colorMap = zeros(3,64);
    colorMap(3,:) = linspace(0,1,64);
    colorMap = colorMap';
    image(ax,imgB);
    title(ax,'B channel');
    colormap(colorMap);
    c=colorbar;
    title(c,'\epsilon_y_y');
    c.Ticks = [0,0.33,0.66,1];
    set(c,'TickLabels',{'-0.07','0','0.07','0.14'});
    set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
    axis equal; axis off;
    
    imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_c',num2str(iCluster),'_B','.tif']);
    print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    close(f);
    
    % (4) RGB channel
    f = figure;
    ax = axes;
    image(ax,imgLocal);
    title(ax,'RGB channels');
    set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
    axis equal; axis off;
    
    imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_c',num2str(iCluster),'_rgb','.tif']);
    print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    close(f);
end

% [B] RGB for all clusters
exx_R_Local = uint8(  img_exx_local.*(img_local_cNum>0) *255);
exy_G_Local = uint8(  img_exy_local.*(img_local_cNum>0) *255);
eyy_B_Local = uint8(  img_eyy_local.*(img_local_cNum>0) *255);

bkgd = uint8(zeros(size(exx_R_Local)));
imgR = cat(3,exx_R_Local, bkgd, bkgd);
imgG = cat(3,bkgd, exy_G_Local, bkgd);
imgB = cat(3,bkgd, bkgd, eyy_B_Local);
% imgLocal = cat(3,imgLocal,imgLocal,imgLocal);
imgLocal = cat(3,exx_R_Local,exy_G_Local,eyy_B_Local);

% (1) R channel
f = figure;
ax = axes;
colorMap = zeros(3,64);
colorMap(1,:) = linspace(0,1,64);
colorMap = colorMap';
image(ax,imgR);
title(ax,'R channel');
colormap(colorMap);
c=colorbar;
title(c,'\epsilon_x_x');
c.Ticks = [0,0.33,0.66,1];
set(c,'TickLabels',{'-0.14','-0.07','0','0.07'});
set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
axis equal; axis off;

imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_gR','.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);

% (2) G channel
f = figure;
ax = axes;
colorMap = zeros(3,64);
colorMap(2,:) = linspace(0,1,64);
colorMap = colorMap';
image(ax,imgG);
title(ax,'G channel');
colormap(colorMap);
c=colorbar;
title(c,'\epsilon_x_y');
c.Ticks = [0,0.5,1];
set(c,'TickLabels',{'-0.07','0','0.07'});
set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
axis equal; axis off;

imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_gG','.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);

% (3) B channel
f = figure;
ax = axes;
colorMap = zeros(3,64);
colorMap(3,:) = linspace(0,1,64);
colorMap = colorMap';
image(ax,imgB);
title(ax,'B channel');
colormap(colorMap);
c=colorbar;
title(c,'\epsilon_y_y');
c.Ticks = [0,0.33,0.66,1];
set(c,'TickLabels',{'-0.07','0','0.07','0.14'});
set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
axis equal; axis off;

imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_gB','.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);

% RGB channel
f = figure;
ax = axes;
image(ax,imgLocal);
title(ax,'RGB channels');
set(ax,'fontsize',18,'xTickLabel',{''},'yTickLabel',{''});
axis equal; axis off;

imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_grgb','.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
close(f);












%% [2] the cluster number map of the grain of interest at different strain levels
close all;
ID_local = ID(indR_min:indR_max, indC_min:indC_max);
colorMap = parula(3);   % This may need to be adjusted to 5 in some cases
colorMap = lines(4);
colorMap(2,:) = [];
colorMap = [.98 .98 .98; colorMap];

for ii = 1:length(iE_list)
    
    iE = iE_list(ii);
    iC = iC_list(ii);
    
    cMapA = cNumMaps{iE}(indR_min:indR_max, indC_min:indC_max);
    cMapA_cleaned = cNumMaps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
    
    cMapA(ID_local~=ID_current) = 0;  % cluster number just this grain
    cMapA_cleaned(ID_local~=ID_current) = 0;
    cNums = unique(cMapA(:));
    cNums(isnan(cNums)) = [];
    n = length(cNums);
    
    % cluster number map, not cleaned
    [f,a,c] = myplot(cMapA);
    colormap(colorMap(1:n,:));
    caxis(caxis+[-0.5 0.5]);
    c.Ticks = 1:n;
    set(c,'limits',[0.5,n-0.5]);
    title(a,'Cluster Label Map');
    set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
    axis equal;
    imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_clusterNumMap.tif']);
    print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    close(f);
    
    % cluster number map, not cleaned
    [f,a,c] = myplot(cMapA_cleaned);
    colormap(colorMap(1:n,:));
    caxis(caxis+[-0.5 0.5]);
    c.Ticks = 1:n;
    set(c,'limits',[0.5,n-0.5]);
    title(a,'Cluster Label Map');
    set(a,'fontsize',24,'xTickLabel',{''},'yTickLabel',{''});
    axis equal;
    imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_clusterNumMap_cleaned.tif']);
    print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    close(f);
end


%% [2 text summary] summary of overlap at different strain levels, overlap of all clusters.  [No plot, just display some value.]
clc
for iE = iE_list(1:end-1)
    
    struA = struCell{iE};   % pre
    struP = struCell{iE+1}; % post
    
    iS = find(arrayfun(@(x) x.gID == ID_target,stru));  % for debugging
    
    ID_current = gIDwithTrace(iS);
    
    ind_local = ismember(ID, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
    indC_min = find(sum(ind_local, 1), 1, 'first');
    indC_max = find(sum(ind_local, 1), 1, 'last');
    indR_min = find(sum(ind_local, 2), 1, 'first');
    indR_max = find(sum(ind_local, 2), 1, 'last');
    
    ID_local = ID(indR_min:indR_max, indC_min:indC_max);
    
    cMapA = cNumMaps{iE}(indR_min:indR_max, indC_min:indC_max);
    cMapP = cNumMaps{iE+1}(indR_min:indR_max, indC_min:indC_max);
    cMapA_cleaned = cNumMaps_cleaned{iE}(indR_min:indR_max, indC_min:indC_max);
    cMapP_cleaned = cNumMaps_cleaned{iE+1}(indR_min:indR_max, indC_min:indC_max);
    
    cMapA(ID_local~=ID_current) = 0;  % cluster number just this grain
    cMapP(ID_local~=ID_current) = 0;
    cMapA_cleaned(ID_local~=ID_current) = 0;
    cMapP_cleaned(ID_local~=ID_current) = 0;
    
    % match clusters, and find the one with area grown.
    cOverlap_noClean = [];
    cOverlap = [];
    for ii=1:length(struA(iS).cLabel)
        cNumA = struA(iS).cLabel(ii);
        for jj=1:length(struP(iS).cLabel)
            cNumP = struP(iS).cLabel(jj);
            
            cOverlap_noClean(ii,jj) =sum(sum((cMapA==cNumA)&(cMapP==cNumP)));
            cOverlap(ii,jj) =sum(sum((cMapA_cleaned==cNumA)&(cMapP_cleaned==cNumP)));   % use cleaned map to calculate how good they overlap !!! But NOT the actual cluster size ------------------
            
        end
    end
    
    volA = sum(cOverlap,2);
    volA_noClean = sum(cOverlap_noClean,2);
    
    % this is more precise:
    volA_noClean = struA(iS).cVol;
    volA = struA(iS).cVolCleaned;
    
    % vrFwd = volP./volA;
    overlapPctA_noClean = cOverlap_noClean./volA_noClean;
    overlapPctA = cOverlap./volA;
    
    % do an additional clean up of cOverlap, ignore those which only have <15% overlap with post-cluster
    % Here, can disable for the purpose of the paper !!! don't want to show 0.  Instead, show actual overlap !!!   -------------- !!!! ----- !!!!! -----
    %     cOverlap(overlapPctA<0.15) = 0;
    
    [cFrom,cTo] = hungarian_assign(max(cOverlap_noClean(:))-cOverlap_noClean);
    link_noClean = false(size(cOverlap_noClean));
    for ii = 1:length(cFrom)
        if (cFrom(ii)>0)&&(cTo(ii)>0)
            link_noClean(cFrom(ii),cTo(ii)) = true;
        end
    end
    
    [cFrom,cTo] = hungarian_assign(max(cOverlap(:))-cOverlap);
    link = false(size(cOverlap));
    for ii = 1:length(cFrom)
        if (cFrom(ii)>0)&&(cTo(ii)>0)
            link(cFrom(ii),cTo(ii)) = true;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    disp('----------------------')
    disp(['strain: ', num2str(iE), ', grain: ',num2str(ID_current)]);
    disp('Non cleaned data: ');
    disp('cluster overlap non-cleaned, cOverlap_noClean: ')
    disp(cOverlap_noClean)
    disp('link noClean: ')
    disp(link_noClean)
    disp('pct overlaped in cluster at early strain, overlapPctA_noClean: ')
    disp(overlapPctA_noClean)
    disp('cluster size, volA_noClean')
    disp(volA_noClean)
    
    disp('cluster overlap, cOverlap: ')
    disp(cOverlap)
    disp('link noClean: ')
    disp(link)
    disp('pct overlaped in cluster at early strain, overlapPctA: ')
    disp(overlapPctA)
    disp('cluster size, volA')
    disp(volA)
    
    if(iE == iE_list(end-1))
        
        disp('----------------------')
        disp(['strain: ', num2str(iE+1), ', grain: ',num2str(ID_current)]);
        
        disp('cluster size, volA_noClean')
        disp(struP(iS).cVol)
        
        disp('cluster size, volA')
        disp(struP(iS).cVolCleaned)
    end
    
end

%% [2 something else] Plot cluster size/volume evolution
close all;
stru = struCell{iE_target};
volEvo = stru(iS).volEvo(iC_target,iE_list);
volEvoCleaned = stru(iS).volEvoCleaned(iC_target,iE_list);

f = figure;
ax = axes;
plot(iE_list,volEvo,'-bo','linewidth',1.5,'markersize',8);
xlabel('Strain Level');
ylabel('Cluster Size, pixels');
title('Cluster Size Evolution');
set(ax,'fontsize',18,'xlim',[1.9, 5.1],'ylim',[0,65000]);
annotation(f,'textbox', [0.15 0.25 0.5 0.075], 'String',{'Cluster # 1 at Strain Level 2'}, 'FontSize',16,'LineStyle','none');
box on;
imgName = (['g',num2str(ID_current),'_s',num2str(iE_target),'_c',num2str(iC_target),'_sizeEvolution.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder

f = figure;
ax = axes;
hold on;
plot(iE_list,volEvo,'-bo','linewidth',1.5,'markersize',8);
plot(iE_list,volEvoCleaned,'-md','linewidth',1.5,'markersize',8);
xlabel('Strain Level');
ylabel('Cluster Size, pixels');
title('Cluster Size Evolution');
set(ax,'fontsize',18,'xlim',[1.9, 5.1],'ylim',[0,65000]);
legend({'uncleand','cleaned'},'Location','best');
annotation(f,'textbox', [0.15 0.25 0.5 0.075], 'String',{'Cluster # 1 at Strain Level 2'}, 'FontSize',16,'LineStyle','none');
box on;
imgName = (['g',num2str(ID_current),'_s',num2str(iE_target),'_c',num2str(iC_target),'_sizeEvolution_withClean.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder





%% [3] Evaluate the applied cluster result, by 3D scatter plot
% close all;
colors = lines(7);
colors(2,:) = [];

figure; 
a = axes;
hold on;
plot3(a,stru(iS).tStrain(:,1),stru(iS).tStrain(:,2),stru(iS).tStrain(:,3),'.k','markersize',24);

step = 20;
for ii = 1:3
    xx = exx_local(clusterNumMapLocal==ii);
    yy = exy_local(clusterNumMapLocal==ii);
    zz = eyy_local(clusterNumMapLocal==ii);
%     plot3(xx(1:step:end),yy(1:step:end),zz(1:step:end),'.','color',colors(ii,:));
    scatter3(a,xx(1:step:end),yy(1:step:end),zz(1:step:end),10,'MarkerEdgeColor','none','MarkerFaceColor',colors(1,:),'MarkerFaceAlpha',0.33);
end
xlabel('\epsilon_x_x');
ylabel('\epsilon_x_y');
zlabel('\epsilon_y_y');
set(a,'fontsize',18,'xlim',[-0.21,0.14],'ylim',[-0.07,0.14],'zlim',[-0.14 0.14])

view([0 0 1]);
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_distribution_xy.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder

view([0 1 0]);
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_distribution_xz.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder

view([1 0 0]);
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_distribution_yz.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder

view(25,45);
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_distribution.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder

%% Illustrate the [process] of how to determine [nCluster] to use.  Here are some comparisons of using different nCluster to try. 
close all;
colors = lines(7);
colors(2,:) = [];
rng(0);     % adjust this if you don't like the cluster color

maxCluster = 5;

ind = find((ID_local==ID_current)); % ind = find((ID_local==ID_current)&(~isnan(exx_local))&(~isnan(exy_local))&(~isnan(eyy_local)));
exx_t = exx_local(ind);
exy_t = exy_local(ind);
eyy_t = eyy_local(ind);
data_t = [exx_t(:), exy_t(:), eyy_t(:)];

% first, predict centroid
[~,ind_centroid_initial] = max(stru(iS).tSF);
centroid_initial = stru(iS).tStrain(ind_centroid_initial,:);

% ======================= kmeans, determine optimum number of clusters ====================================
% % The way it was actually used.
% nPoints = 8100;
% ind_reduce = ~isnan(sum(data_t,2));
% data_reduce = data_t(ind_reduce,:);
% reduce_ratio = ceil(size(data_reduce,1)/nPoints);
% data_reduce = data_reduce(1:reduce_ratio:end,:);

% To fit here:
data_reduce = data_t(1:step:end,:);

if(~isempty(data_reduce))
    % compare the silhouette, by actually do kmeans on down-sampled samples.
    disp(['ID=',num2str(ID_current)]);
    clear wssd  score_avg  score_cluster_mean  score_cluster_neg_sum  score_cluster_mean_min  score_neg_sum;
    %         score_min = -1*ones(1, maxCluster);
    score_neg_sum = -inf*ones(1, maxCluster);
    nRep = 1;
    c0 = kmeans_pp_init(data_reduce,maxCluster,nRep,centroid_initial);
    for nc = 2:maxCluster
        [idx, centroid, sumd] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',1000,'start',c0(1:nc,:,:));   % 'correlation' distance not good.
        sil_score = silhouette(data_reduce,idx);
        
        wssd(nc) = mean(sumd);
        score_avg(nc) = nanmean(sil_score); % avg score for the condition of nc clusters
        
        for ii=1:nc
            sil_this_cluster = sil_score(idx==ii);
            score_cluster_mean{nc}(ii) = mean(sil_this_cluster); % silhouette for each cluster
            score_cluster_neg_sum{nc}(ii) = sum(sil_this_cluster(sil_this_cluster<0));
        end
        
        score_cluster_mean_min(nc) = min(score_cluster_mean{nc});
        score_neg_sum(nc) = sum(score_cluster_neg_sum{nc});
        
        % the Silhouette figures:
        % (1) just the figure
        figure; silhouette(data_reduce,idx);
        set(gca,'fontsize',18,'xlim',[-0.5,1]);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_silhouette_nc_',num2str(nc),'.tif']);
        print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder

        % (2) with average sihouette
        figure; silhouette(data_reduce,idx);
        set(gca,'fontsize',18,'xlim',[-0.5,1]);
        hold on; fplot(@(x) (x-score_avg(nc))*2^64,'--r','linewidth',2);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_silhouette_nc_',num2str(nc),'_withAvg.tif']);
        print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
        
        % (3) scatter plot show the distribution of the clusters
        figure;
        a = axes;
        hold on;
        plot3(a,stru(iS).tStrain(:,1),stru(iS).tStrain(:,2),stru(iS).tStrain(:,3),'.k','markersize',24);
        
         for ii = 1:nc
            xx = data_reduce(idx==ii,1);
            yy = data_reduce(idx==ii,2);
            zz = data_reduce(idx==ii,3);
            %     plot3(xx(1:step:end),yy(1:step:end),zz(1:step:end),'.','color',colors(ii,:));
            scatter3(a,xx,yy,zz,10,'MarkerEdgeColor','none','MarkerFaceColor',colors(ii,:),'MarkerFaceAlpha',0.33);
        end
        xlabel('\epsilon_x_x');
        ylabel('\epsilon_x_y');
        zlabel('\epsilon_y_y');
        set(a,'fontsize',18,'xlim',[-0.2,0.05],'ylim',[-0.07,0.09],'zlim',[-0.09 0.14]);
        
        view([0 0 1]);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_scatter_nc_',num2str(nc),'_xy.tif']);
        print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
        
        view([0 1 0]);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_scatter_nc_',num2str(nc),'_xz.tif']);
        print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
        
        view([1 0 0]);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_scatter_nc_',num2str(nc),'_yz.tif']);
        print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
        
        view(25,45);
        imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_scatter_nc_',num2str(nc),'.tif']);
        print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    end
    %         [~,nCluster] = max(score_cluster_mean_min);
    [~,nCluster] = max(score_neg_sum);
    disp([char(9), 'ID= ',num2str(ID_current),', nCluster=',num2str(nCluster)]);
    
    sNC.gID = ID_current;
    sNC.nCluster = nCluster;
    sNC.wssd = wssd;
    sNC.score_avg = score_avg;
    sNC.score_cluster_mean = score_cluster_mean;
    sNC.score_cluster_neg_sum = score_cluster_neg_sum;
    sNC.score_cluster_mean_min = score_cluster_mean_min;
    sNC.score_neg_sum = score_neg_sum;
end

%% [caution!] when satisfied, use this template to save image
% view(25,45);
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_eij_scatter_nc_3.tif']);
imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_silhouette_nc_3.tif']);
print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
    

%% continue on 2018-03-21, examine pointwise distance to twin strain
% try using K=2-5 clusters, plot
% (1) pointwise distance of strain components to highest-SFed twin center.
% This could be used to make a contour plot of equi-distance, but it turned
% out hard to see, so just use this pointwise distance map.
% (2) The histogram of this distance. The boundary of twin vs non-twin can
% be guessed from the position between peaks. -- Important thing is, the
% boundary position does not seem to be constant.

% First, iE of interest ----------------------------- !!! -------------------------------------
iE = 5;
fname = [f1,'_all_grain_',num2str(ID_current),'_local_map.mat'];
load(fullfile(grainDataPath,fname),'data');
s = data(iE);

indR_min = s.indR_min;
indR_max = s.indR_max;
indC_min = s.indC_min;
indC_max = s.indC_max;
boundaryTF_local = s.boundaryTF_local;
x_local = s.x_local;
y_local = s.y_local;
ID_current = s.ID_current;
ID_local = s.ID_local;
exx_local = s.exx_local;
exy_local = s.exy_local;
eyy_local = s.eyy_local;
sigma_local = s.sigma_local;

% summarize the distribution of exx data at different iE level
% figure;histogram(exx_local(:));
% figure;histogram(exy_local(:));
% figure;histogram(eyy_local(:));
edges_xx = linspace(-0.2,0.05,100);
edges_xy = linspace(-0.04,0.06,100);
edges_yy = linspace(-0.04,0.06,100);
exx_distribution{iE} = histcounts(exx_local(:),edges_xx);
exy_distribution{iE} = histcounts(exy_local(:),edges_xx);
eyy_distribution{iE} = histcounts(eyy_local(:),edges_xx);


% close all;
colors = lines(7);
colors(2,:) = [];
rng(0);     % adjust this if you don't like the cluster color

maxCluster = 1;

ind = find((ID_local==ID_current)); % ind = find((ID_local==ID_current)&(~isnan(exx_local))&(~isnan(exy_local))&(~isnan(eyy_local)));
exx_t = exx_local(ind);
exy_t = exy_local(ind);
eyy_t = eyy_local(ind);
data_t = [exx_t(:), exy_t(:), eyy_t(:)];

if 0 % re-cluster and analyze.  But turned out not able to make contour plot well.  So just keep code and disable.
    % first, predict centroid
    [~,ind_centroid_initial] = max(stru(iS).tSF);
    centroid_initial = stru(iS).tStrain(ind_centroid_initial,:);
    
    % To fit here:
    data_reduce = data_t(1:end,:);
    
    if(~isempty(data_reduce))
        % compare the silhouette, by actually do kmeans on down-sampled samples.
        disp(['ID=',num2str(ID_current)]);
        nRep = 1;
        c0 = kmeans_pp_init(data_reduce,maxCluster,nRep,centroid_initial);
        
        pointWiseDistMapLocal = zeros(size(exx_local));      % record raw clusterNumberMap
        pointWiseDist = pdist2(centroid_initial, data_reduce);
        pointWiseDistMapLocal(ind) = pointWiseDist;
        
        myplot(pointWiseDistMapLocal); caxis([0 0.05]);
        figure;histogram(pointWiseDist);
        
        for nc = 2:maxCluster
            [idx, centroid, sumd] = kmeans(data_reduce, nc, 'Distance','sqeuclidean','MaxIter',1000,'start',c0(1:nc,:,:));   % 'correlation' distance not good.
            clusterNumMapLocal = zeros(size(exx_local));      % record raw clusterNumberMap
            clusterNumMapLocal(ind) = idx;
            
            %         imgName = (['s',num2str(iE),'_g',num2str(ID_current),'_silhouette_nc_',num2str(nc),'_withAvg.tif']);
            %         print(fullfile(grainDataPath,imgName),'-dtiff');   % to parent folder
            
        end
    end
    
end
disp(['strain level: ',num2str(iE),' analyzed']);






