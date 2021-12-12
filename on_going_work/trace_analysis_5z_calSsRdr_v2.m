% Zhe Chen, 2015-12-6
%
% For each grain, for each slip system, calculate its theoretical RDR.
% Add a field 'ssDispRatio' into the structure.
%
% chenzhe, 2017-06-07.
% changed the name 'ssDispRatio' to 'ssRDR'.
% use the new function trace_analysis_TiMgAl().  But also calculated ssRDR
% in step-2.  So this function is actually no longer needed.

[traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
traceStruct = load([traceStructPath,'\',traceStructFile]);
traceStruct = traceStruct.traceStruct;

previousData = 'Ti7Al_#B6_EbsdToSemForTraceAnalysis';   % The data, that you saved, that contains ebsdToSEM data.
load(previousData);
sampleName=sampleName; sampleMaterial=sampleMaterial; stressTensor=stressTensor; % These are some of the prevously loaded settings.

for iS = 1:length(traceStruct)
    idThis = traceStruct(iS).ID;
    
    ind_current = find(idThis == gID);    % an index of row
    phi1_current = gPhi1(ind_current);
    phi_current = gPhi(ind_current);
    phi2_current = gPhi2(ind_current);
%     [sf_mat, ~, burgersXY] =
%     trace_analysis_Ti(phi1_current,phi_current,phi2_current,-90,180,0,0,0,0)    % choose 'Ti' or 'Mg' ---------------  ---------------  ------------  
    [sf_mat, sf_mat_sorted, burgersXY] = trace_analysis_TiMgAl([phi1_current,phi_current,phi2_current],[-90,180,0],[0,0,0],stressTensor,sampleMaterial,'noTwin');
    traceStruct(iS).ssRDR = burgersXY(:,3)';      % changed name 'ssDispRatio' to 'ssRDR'
    
end

%%
answer = questdlg('save structure?');
if strcmpi(answer, 'Yes')
    try
        save([traceStructPath,'\traceStruct_withSsRDR_rename_it'],'traceStruct','-append');
    catch
        save([traceStructPath,'\traceStruct_withSsRDR_rename_it'],'traceStruct');
    end
end

