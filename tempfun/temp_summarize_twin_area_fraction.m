close all;

iE=6;
strainFile = [dicPath,'\',f2,STOP{iE+B}]; disp(strainFile);
load(strainFile,'exx','exy','eyy');
name_result_modified = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_modified.mat'];
load([saveDataPath,name_result_modified],'twinMap','sfMap','disSimiMap','twinMap_2','shearMap','sfMap_2','costMap','clusterNumMap','stru');

myplot(X,Y,exx,boundaryTFB);
% myplot(X,Y,clusterNumMap,boundaryTFB);

% myplot(X,Y,disSimiMap,boundaryTFB);
myplot(X,Y,twinMap,boundaryTFB); caxis([18,24]);
% myplot(X,Y,sfMap,boundaryTFB);

% myplot(X,Y,shearMap, boundaryTFB);
% myplot(X,Y,abs(shearMap-0.1289), boundaryTFB); title('shear diff');
myplot(X,Y,twinMap_2, boundaryTFB); caxis([18,24]);
% myplot(X,Y,sfMap_2, boundaryTFB);
% myplot(X,Y,costMap, boundaryTFB);