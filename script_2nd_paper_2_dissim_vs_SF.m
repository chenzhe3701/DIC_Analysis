

% just look at dissimilarity vs SF.
% chenzhe, 2018-02-09

clear;
addChenFunction;

% looks like have to include this part to read the sample name.
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');

saveFig = 1;
if saveFig
    saveFigurePath = [uigetdir('C:\Users\ZheChen\Desktop','choose a path to save the figures'),'\'];
end

% load previous data and settings
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);

gIDwithTrace = gID(~isnan(gExx));

% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------

STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

%% plot dissimilarity vs schmid factor 
postAnalysis = 0;
SF = [];
Dis = [];
Twinned = [];
for iE = iE_start:iE_stop
    
    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];    % with or without 'fly_rss'
    if postAnalysis
        name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];       % use this after labeled
    end
    load([saveDataPath,name_result_on_the_fly]);
    
    rtss_cCen_tStrain = [];
    rtss_cMed_tStrain = [];
    schmidfactor = [];
    trueTwin = [];
    for iS=1:length(stru)
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
            if postAnalysis
                % if it is later identified as a twin
                if stru(iS).trueTwin(iCluster)>0
                    trueTwin = [trueTwin,1];
                else
                    trueTwin = [trueTwin,2];
                end
                rtss_cCen_tStrain  = [rtss_cCen_tStrain, stru(iS).dis(iCluster)];
                schmidfactor = [schmidfactor, stru(iS).sf(iCluster)];
            else
                allAvailable = 0;   
                if allAvailable
                    rtss_cCen_tStrain  = [rtss_cCen_tStrain, pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain)];
                    % rtss_cMed_tStrain  = [rtss_cMed_tStrain, pdist2(stru(iS).cMed(iCluster,:), stru(iS).tStrain)];
                    schmidfactor = [schmidfactor, stru(iS).tSF];
                else
                    d = pdist2(stru(iS).cCen(iCluster,:), stru(iS).tStrain);
                    [d,ind] = min(d);
                    rtss_cCen_tStrain  = [rtss_cCen_tStrain,d];
                    % rtss_cMed_tStrain  = [rtss_cMed_tStrain, pdist2(stru(iS).cMed(iCluster,:), stru(iS).tStrain)];
                    schmidfactor = [schmidfactor, stru(iS).tSF(ind)];
                end
            end
        end
    end
    Dis{iE} = rtss_cCen_tStrain;
    SF{iE} = schmidfactor;
    Twinned{iE} = trueTwin;
end


%% plot
close all;
titleString = {'','Global Uniaxial Strain = -0.6%','Global Uniaxial Strain = -1.2%',...
    'Global Uniaxial Strain = -2.1%', 'Global Uniaxial Strain = -3.7%'};
titleString = {'',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.004',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.012',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.023',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.039'};

for iE = iE_start:iE_stop
    
    if postAnalysis
       figure; hold on;
       plot(Dis{iE}(Twinned{iE}==2),SF{iE}(Twinned{iE}==2),'.b','markersize',8);
       plot(Dis{iE}(Twinned{iE}==1),SF{iE}(Twinned{iE}==1),'.r','markersize',8);
       legend({'Non-twinned','Twinned'},'Position',[0.3 0.25 0.3 0.14]);
    else
        figure;plot(Dis{iE},SF{iE},'.','markersize',8);
    end
    
    title(titleString{iE},'fontweight','normal','fontsize',18);
    xlabel('Dissimilarity');
    xlabel('\psi^D_{min}');
    ylabel('Schmid Factor'); 
    ylabel('Schmid Factor (m)');
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.8]);
    
    imgName = (['s',num2str(iE),'_Dis_vs_SF.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
        
        imgNameA = (['s',num2str(iE),'_Dis_vs_SF_a.tif']);
        title(['Strain Level: ',num2str(STOP{iE+B})],'fontweight','normal','fontsize',18);
        print(fullfile(saveFigurePath,imgNameA),'-dtiff');   % to parent folder
        close all;
    end
    
    % zoom-in view
    asp = daspect;
    if postAnalysis
       figure; hold on;
       plot(Dis{iE}(Twinned{iE}==2),SF{iE}(Twinned{iE}==2),'.b','markersize',10);
       plot(Dis{iE}(Twinned{iE}==1),SF{iE}(Twinned{iE}==1),'.r','markersize',10);
    else
        figure;plot(Dis{iE},SF{iE},'.','markersize',10);
    end
    hold on; fplot(@(x) 7*x+0.15,[0 0.05],'--k','linewidth',5);
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.08],'ylim',[0.2, 0.5]);
    daspect(asp);
    
    
    imgName = (['s',num2str(iE),'_Dis_vs_SF_zoom.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
        close all;
    end
%     figure;plot(rtss_cMed_tStrain,schmidfactor,'.');
%     title(['Strain level: ',num2str(STOP{iE+B})]);
%     xlabel('RootSS Cluster Medium - Twin Strain');
%     ylabel('Schmid Factor');
%     
%     figure; histogram(rtss_cCen_tStrain - rtss_cMed_tStrain);
    
end