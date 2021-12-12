% clear;
% addChenFunction;
% dicPath = uigetdir('D:\WE43_T6_C1_insitu_compression\stitched_DIC','pick DIC directory, which contains the stitched DIC data for each stop');
% dicFiles = dir([dicPath,'\*.mat']);
% dicFiles = struct2cell(dicFiles);
% dicFiles = dicFiles(1,:)';
% 
% % looks like have to include this part to read the sample name.
% [fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
% load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');
% 
% % load previous data and settings
% saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
% saveDataPathInput = saveDataPath;
% load([saveDataPath,sampleName,'_traceAnalysis_WS_settings.mat']);
% if ~strcmpi(saveDataPath,saveDataPathInput)
%     disp('Input saveDataPath is different from that saved in setting.mat file. Check files or code.');
%     return;
% end
% try
%     load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis']);
% catch
%     load([saveDataPath,sampleName,'_EbsdToSemForTraceAnalysis_GbAdjusted']);
% end
% 
% gIDwithTrace = gID(~isnan(gExx));
% 
% % modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
% 
% STOP = {'0','1','2','3','4','5','6','7'};
% B=1;    % 0-based B=1.  1-based B=0.
% iE_start = 2;   % elongation levels to analyze. 0-based.
% iE_stop = 6;
% 
% % file name prefixes
% f1 = 'WE43_T6_C1_s';
% f2 = '_';
% 
% neighbor_elim = 1;          % don't consider this ID as neighbor. For example, ID = 1 or 0 means bad region.
% twinTF_text = 'twin';        % do you want to analyze twin? Use things like 'twin' or 'notwin'
% 
%%
clear;

load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\WE43_T6_C1_EbsdToSemForTraceAnalysis.mat','ID','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2');
save('WE43_T6_C1_EBSD_organized.mat','ID','gID','gPhi1','gPhi','gPhi2','phi1','phi','phi2');

load('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\WE43_T6_C1_EbsdToSemForTraceAnalysis_GbAdjusted.mat','ID','gID','gPhi1','gPhi','gPhi2');
save('WE43_T6_C1_EBSD_GbAdjusted_organized.mat','ID','gID','gPhi1','gPhi','gPhi2');



