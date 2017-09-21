% ZheChen 2015-12-3
% identify Active SS based on pre-defined criterion, for example
% 
% (1) angle to observed direction is < 6 degree.
% (2) among those, angle_diff_normalized / weighed_SF, choose the smallest one.
% 
% Using the above criterion, if only 1 ss matches well, then the label 'TC' (for 'traceCertain') is set to be '1'
%
% If ONLY 1 ss lie within (6) degrees of labeled trace, then the label 'TC1' is set to be '1'
%
% in summary, in the current version, (two) fields labeling ss
% identification quality is created.
%
% chenzhe 2017-06-07.
% Add Al.
%
% new fields after this step: 'ssActivated','TC','TC1'

[traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
traceStruct = load([traceStructPath,'\',traceStructFile]);
traceStruct = traceStruct.traceStruct;

previousData = 'T5_#7_EbsdToSemForTraceAnalysis';   % The data, that you saved, that contains ebsdToSEM data.
load(previousData);
sampleName=sampleName; sampleMaterial=sampleMaterial; stressTensor=stressTensor; % These are some of the prevously loaded settings.

%%
% ssWeights can be modified later -------------------------------------------------- 
if any(strcmpi(sampleMaterial,{'Mg','Magnesium'}))
    ssWeight = [1 1 1, 5 5 5, 5 5 5 5 5 5, 10 10 10 10 10 10, 2 2 2 2 2 2];    
elseif any(strcmpi(sampleMaterial,{'Ti','Titanium'}))
    ssWeight = [1.25 1.25 1.25, 1 1 1, 7.4 7.4 7.4 7.4 7.4 7.4, 24 24 24 24 24 24 24 24 24 24 24 24, 1.25 1.25 1.25 1.25 1.25 1.25];
elseif any(strcmpi(sampleMaterial,{'Al','Aluminum'}))
    ssWeight = [1 1 1 1 1 1 1 1 1 1 1 1]; % Al without twinning
end

for iS = 1:length(traceStruct)   % iRow, for each grain
    traceStruct(iS).ssActivated = zeros(size(traceStruct(iS).ssActivated));  % Note: clear for reanalyze !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for iTrace = 1:traceStruct(iS).nTraces   % iTrace, for each observed trace
        
        angle_diff = ceil(abs(traceStruct(iS).ssTraceDir - traceStruct(iS).traceDir(iTrace)));    % ss theoretical trace direction - observed direction
        minAngleDiffSatisfied = find(angle_diff <= 6);    % has to be smaller than (some, 6 for now) degree difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        notActiveYet = find(traceStruct(iS).ssActivated==0);     % has not been assigned
        
        candidateSS = intersect(minAngleDiffSatisfied, notActiveYet);
        
        % TC: traceCertain: if one ss match twice better than 2nd choice
        % TC1: if only one lie within 6 degrees.
        if ~isempty(candidateSS)
            % Both initialize as zero.
            traceStruct(iS).TC(iTrace) = 0;
            traceStruct(iS).TC1(iTrace) = 0;
        
            angle_diff = angle_diff(candidateSS);
            angle_diff_ratio = angle_diff/min(angle_diff);     % angle_diff normalize by the smallest one
            if sum(angle_diff_ratio <5) < 2      % if one trace matches very well than all other, --------------------- this affects the result.  If all not clear, maybe change this larger
                [~, ind_temp] = min(angle_diff_ratio);
                activatedSS = candidateSS(ind_temp);
                traceStruct(iS).TC(iTrace) = 1;      % it is certain that it activated (current criterion is only based on direction comparison)
                if length(candidateSS) == 1
                    traceStruct(iS).TC1(iTrace) = 1;
                end
            else
                vote = angle_diff_ratio ./ ( (traceStruct(iS).sf(candidateSS)) ./ ssWeight(candidateSS) );
                [~, ind_temp] = min(vote);
                activatedSS = candidateSS(ind_temp);
            end
            traceStruct(iS).ssActivated(activatedSS)  = iTrace;     % if activated, label it with corresponding trace#
        end
    end
end
%%
answer = questdlg('save structure?');
if strcmpi(answer, 'Yes')
    try
        save('traceStruct_ssIdentified_rename_it','traceStruct','-append');
        disp('save append');
    catch
        save('traceStruct_ssIdentified_rename_it','traceStruct');
    end
end

% save('traceStruct_new','traceStruct');

