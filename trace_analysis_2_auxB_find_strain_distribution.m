
% find how experimental strain values distribute, to get an idea how to
% smooth data
%
% chenzhe, 2018-02-04

dicPath = uigetdir('','pick DIC directory, which contains the stitched DIC data for each stop');
dicFiles = dir([dicPath,'\*.mat']);
dicFiles = struct2cell(dicFiles);
dicFiles = dicFiles(1,:)';
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
f_prefix= '_';
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

qt_exx=[];
qt_exy=[];
qt_eyy=[];

qts = cdf('normal',[-5:0.5:5], 0, 1);   % x=multiples of sigma, 0=mean, 1=sigma 
%%
for iE=2
strainFile = [dicPath,'\',f_prefix,STOP{iE+B}]; disp(strainFile);
load(strainFile,'exx','exy','eyy','sigma');     
    exx(sigma==-1)=nan;
    exy(sigma==-1)=nan;
    eyy(sigma==-1)=nan;
    qt_exx = [qt_exx; quantile(exx(:),qts)];
    qt_exy = [qt_exy; quantile(exy(:),qts)];
    qt_eyy = [qt_eyy; quantile(eyy(:),qts)];
    
    exx_mean = nanmean(exx(:));
    exx_std = nanstd(exx(:));
    exy_mean = nanmean(exy(:));
    exy_std = nanstd(exy(:));
    eyy_mean = nanmean(eyy(:));
    eyy_std = nanstd(eyy(:));
    

end

disp(qt_exx);
disp(qt_exy);
disp(qt_eyy);
figure;histogram(exx(:),'binlimits',[-0.4 0.4]);
figure;histogram(exy(:),'binlimits',[-0.4 0.4]);
figure;histogram(eyy(:),'binlimits',[-0.4 0.4]);

% a = histcounts(exx(:),'binedges',[-1:0.01:1]);

% exx((exx<exx_mean-1.5*exx_std)|(exx>exx_mean+1.5*exx_std))=nan;
% exy((exy<exy_mean-1.5*exy_std)|(exy>exy_mean+1.5*exy_std))=nan;
% eyy((eyy<eyy_mean-1.5*eyy_std)|(eyy>eyy_mean+1.5*eyy_std))=nan;
% 
% myplot(exx);
% myplot(exy);
% myplot(eyy);


%%
save([saveDataPath,'strain_quantiles_every_p5.mat'],'qt_exx','qt_exy','qt_eyy');
