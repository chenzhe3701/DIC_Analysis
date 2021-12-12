
function [] = cluster_ID_colorbar(lMin, lMax, nN, bgColor)
% make colorbar for things like cluster-ID-map
% lMin = smallest ID
% lMax = largest ID
% N = number of clusters
% chenzhe, 2019-04-25

delta = (lMax - lMin)/nN;

cMax = lMax + delta;
cMin = lMin - 3*delta;
zMax = lMax + delta;
zMin = lMin - delta;

cMap = parula(nN+1);
if exist('bgColor','var')
    cMap(1,:) = bgColor;
end

colormap(cMap);

caxis([cMin,cMax]);


handle_c = findall(gcf,'type','ColorBar');
if isempty(handle_c)
    display('No ColorBar.');
end

set(handle_c,'limits',[zMin, zMax]);

end