% the confusion map for TP/FP/FN/TN events are summarized in code:
open trace_analysis_3D_clusterToTwin_by_trace_strain.m;

% use the following to choose files to compare
[truthFile, truthPath] = uigetfile('D:\p\m\DIC_Analysis\temp_results\WE43_T6_C1_new_variant_map.mat','select the truth results for twin-grain boundary intersection');
[checkFile, checkPath] = uigetfile('D:\p\m\DIC_Analysis\20191230_1756_twinMaps_K.mat','select the results for twin-grain boundary intersection to check');

% Also need to clean up the twinMapCell (as if it is the final map) first before doing summary 
[toCleanFile, toCleanPath] = uigetfile('D:\p\m\DIC_Analysis\20191230_1756_twinMaps_K.mat','select the results to do clean up');