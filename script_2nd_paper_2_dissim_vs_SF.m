

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
%%
colors = lines(7);
figure;hold on;
for ii=1:size(colors,1)
   scatter(ii,ii,'MarkerFaceColor',colors(ii,:)); 
end
%% (1) initial explore, plot dissimilarity vs schmid factor 
allAvailable = 0;   % plot for all available (6 times) or best matched
SF = [];
Dis = [];
for iE = iE_start:iE_stop
    
    name_result_on_the_fly = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_result_on_the_fly.mat'];    % with or without 'fly_rss'
    load([saveDataPath,name_result_on_the_fly]);
    
    rtss_cCen_tStrain = [];
    rtss_cMed_tStrain = [];
    schmidfactor = [];
    trueTwin = [];
    for iS=1:length(stru)
        nCluster = length(stru(iS).cLabel);
        for iCluster = 1:nCluster
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
    Dis{iE} = rtss_cCen_tStrain;
    SF{iE} = schmidfactor;
end


%% (1.2) plot. Dashed: approximate slope with C=7 and at the same approximate position
close all;
% titleString = {'','Global Uniaxial Strain = -0.6%','Global Uniaxial Strain = -1.2%',...
%     'Global Uniaxial Strain = -2.1%', 'Global Uniaxial Strain = -3.7%'};
titleString = {'',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.004',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.012',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.023',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.039'};

for iE = iE_start:iE_stop
    
    figure;plot(Dis{iE},SF{iE},'.b','markersize',8,'Color',colors(1,:));
    
    title(titleString{iE},'fontweight','normal','fontsize',18);
    xlabel('Dissimilarity');
    xlabel('\psi^D_{min}');
    ylabel('Schmid Factor (m)');
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.8]);
    
    imgName = (['s',num2str(iE),'_Dis_vs_SF_explore.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
%         
%         imgNameA = (['s',num2str(iE),'_Dis_vs_SF_a.tif']);
%         title(['Strain Level: ',num2str(STOP{iE+B})],'fontweight','normal','fontsize',18);
%         print(fullfile(saveFigurePath,imgNameA),'-dtiff');   % to parent folder
%         close all;
    end
    
    % zoom-in view
    asp = daspect;
    figure;plot(Dis{iE},SF{iE},'.b','markersize',10,'Color',colors(1,:));
    hold on; fplot(@(x) 7*x+0.15,[0 0.05],'--k','linewidth',4);
    
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.08],'ylim',[0.2, 0.5]);
    daspect(asp);
    
    
    imgName = (['s',num2str(iE),'_Dis_vs_SF_explore_zoom.tif']);
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


%% (2) do post analysis, to plot dissimilarity vs schmid factor 
postAnalysis = 1;
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

%% (2.2) in post analysis, try to fit model
for iE = 2:5
    disp(['---------------', num2str(iE), '-------------']);
    % [b,dev,stats] = mnrfit([Dis{iE}',SF{iE}'],Twinned{iE}');
    [b,dev,stats] = glmfit([Dis{iE}',SF{iE}'],Twinned{iE}'==1,'binomial','link','logit');   % 'logit' vs 'identity'
    display(b);
    display(b(2)/b(3));    
    coef{iE} = b;    
end

%% (2.3) check. Fit, so that a = xb = b0 + x1b1 + x2b2, and prob = logsig(a) = exp(a)/(1+exp(a))
iE = 2;
dis = Dis{iE}(Twinned{iE}==1)';
sf = SF{iE}(Twinned{iE}==1)';
p1 = coef{iE}(1) + dis * coef{iE}(2) + sf * coef{iE}(3);
% logsig(a) = exp(a)/(1+exp(a))
p1 = logsig(p1);    
% this is the same as  = glmval(coef{iE},[dis,sf],'logit');

dis = Dis{iE}(Twinned{iE}==2)';
sf = SF{iE}(Twinned{iE}==2)';
p2 = coef{iE}(1) + dis * coef{iE}(2) + sf * coef{iE}(3);
p2 = logsig(p2);    % = glmval(coef{iE},[dis,sf],'logit')

figure; hold on;
histogram(p1,'FaceColor','r');
histogram(p2,'FaceColor','b');

%% (2.4) plot. solid line: best with C fitted by regression. Dashed line: best with c=7.
close all;
% titleString = {'','Global Uniaxial Strain = -0.6%','Global Uniaxial Strain = -1.2%',...
%     'Global Uniaxial Strain = -2.1%', 'Global Uniaxial Strain = -3.7%'};
titleString = {'',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.004',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.012',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.023',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.039'};

for iE = iE_start:iE_stop
    
    figure; hold on;
    plot(Dis{iE}(Twinned{iE}==2),SF{iE}(Twinned{iE}==2),'.b','markersize',8,'Color',colors(1,:));
    plot(Dis{iE}(Twinned{iE}==1),SF{iE}(Twinned{iE}==1),'.r','markersize',8,'Color',colors(2,:));
    legend({'Non-twinned','Twinned'},'Position',[0.3 0.25 0.3 0.14]);

    title(titleString{iE},'fontweight','normal','fontsize',18);
    xlabel('Dissimilarity');
    xlabel('\psi^D_{min}');
    ylabel('Schmid Factor (m)');
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.8]);
    
    imgName = (['s',num2str(iE),'_Dis_vs_SF_postRegC7.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
    end
    
    % zoom-in view
    asp = daspect;

    figure; hold on;
    plot(Dis{iE}(Twinned{iE}==2),SF{iE}(Twinned{iE}==2),'.','markersize',10,'Color',colors(1,:));
    plot(Dis{iE}(Twinned{iE}==1),SF{iE}(Twinned{iE}==1),'.','markersize',10,'Color',colors(2,:));
    switch iE
        case 2
            fplot(@(x) 5.03*x + 0.2724, [0 0.05],'-','linewidth',4,'Color',[1 0 0]);
            fplot(@(x) 7*x + 0.2255, [0 0.05],'--k','linewidth',4);
        case 3
            fplot(@(x) 6.26*x + 0.2574, [0 0.05],'-','linewidth',4,'Color',[1 0 0]);
            fplot(@(x) 7*x + 0.2217, [0 0.05],'--k','linewidth',4);
        case 4
            fplot(@(x) 6.24*x + 0.2798, [0 0.05],'-','linewidth',4,'Color',[1 0 0]);
            fplot(@(x) 7*x + 0.2186, [0 0.05],'--k','linewidth',4);
        case 5
            fplot(@(x) 3.59*x + 0.3491, [0 0.05],'-','linewidth',4,'Color',[1 0 0]);
            fplot(@(x) 7*x + 0.2517, [0 0.05],'--k','linewidth',4);
    end

    
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.08],'ylim',[0.2, 0.5]);
    daspect(asp);
    
    
    imgName = (['s',num2str(iE),'_Dis_vs_SF_postRegC7_zoom.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
        close all;
    end
end



%% (3) plot after post analysis. Solid: Both slope C and position determined by regression. Regression link with 'logit'.
for iE = 2:5
    disp(['---------------', num2str(iE), '-------------']);
    % [b,dev,stats] = mnrfit([Dis{iE}',SF{iE}'],Twinned{iE}');
    [b,dev,stats] = glmfit([Dis{iE}',SF{iE}'],Twinned{iE}'==1,'binomial','link','logit');   % 'logit' vs 'identity'
    display(b);
    display(b(2)/b(3));    
    coef{iE} = b;    
end

close all;
titleString = {'',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.004',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.012',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.023',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.039'};

for iE = iE_start:iE_stop
    
        figure; hold on;
        plot(Dis{iE}(Twinned{iE}==2),SF{iE}(Twinned{iE}==2),'.b','markersize',8,'Color',colors(1,:));
        plot(Dis{iE}(Twinned{iE}==1),SF{iE}(Twinned{iE}==1),'.r','markersize',8,'Color',colors(2,:));
        
        fplot(@(x) -coef{iE}(2)/coef{iE}(3)*x -coef{iE}(1)/coef{iE}(3), [0 0.05],'-k','linewidth',2.5);
        fplot(@(x) 7*x + 0.2255, [0 0.05],'--k','linewidth',2.5);
        
        legend({'Non-twinned','Twinned','Regression Boundary','Threshold Boundary'},'Position',[0.3 0.3 0.35 0.14],'fontsize',12);

    
    title(titleString{iE},'fontweight','normal','fontsize',18);
    xlabel('Dissimilarity');
    xlabel('\psi^D_{min}');
    ylabel('Schmid Factor (m)');
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.8],'ylim',[-0.6, 0.6]);
        
    imgName = (['s',num2str(iE),'_Dis_vs_SF_Reg_logit.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
    end
end

%% (4) plot after post analysis. Solid: Both slope C and position determined by regression. Regression link with 'identity'.
for iE = 2:5
    disp(['---------------', num2str(iE), '-------------']);
    % [b,dev,stats] = mnrfit([Dis{iE}',SF{iE}'],Twinned{iE}');
    [b,dev,stats] = glmfit([Dis{iE}',SF{iE}'],Twinned{iE}'==1,'binomial','link','identity');   % 'logit' vs 'identity'
    display(b);
    display(b(2)/b(3));    
    coef{iE} = b;    
end

close all;
titleString = {'',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.004',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.012',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.023',...
    '\epsilon\fontsize{12}^G\fontsize{18}=-0.039'};

for iE = iE_start:iE_stop
    
        figure; hold on;
        plot(Dis{iE}(Twinned{iE}==2),SF{iE}(Twinned{iE}==2),'.b','markersize',8,'Color',colors(1,:));
        plot(Dis{iE}(Twinned{iE}==1),SF{iE}(Twinned{iE}==1),'.r','markersize',8,'Color',colors(2,:));
        
        fplot(@(x) -coef{iE}(2)/coef{iE}(3)*x -coef{iE}(1)/coef{iE}(3), [0 0.05],'-k','linewidth',2.5);
        fplot(@(x) 7*x + 0.2255, [0 0.05],'--k','linewidth',2.5);
        
        legend({'Non-twinned','Twinned','Regression Boundary','Threshold Boundary'},'Position',[0.3 0.3 0.35 0.14],'fontsize',12);

    
    title(titleString{iE},'fontweight','normal','fontsize',18);
    xlabel('Dissimilarity');
    xlabel('\psi^D_{min}');
    ylabel('Schmid Factor (m)');
    set(gca,'fontsize',18);
    set(gca,'xlim',[0 0.8],'ylim',[-0.6, 0.6]);
        
    imgName = (['s',num2str(iE),'_Dis_vs_SF_Reg_identity.tif']);
    if saveFig
        print(fullfile(saveFigurePath,imgName),'-dtiff');   % to parent folder
    end
end

close all;