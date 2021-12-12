function traceStruct = trace_analysis_3a_removeTrace_v2(traceStruct, target_id, target_trace_num)

target_ind = find(cell2mat(arrayfun(@(x) (x.ID==target_id), traceStruct, 'uniformoutput',false)));

if traceStruct(target_ind).nTraces >= target_trace_num
    % These are from adding slip traces
    traceStruct(target_ind).nTraces = traceStruct(target_ind).nTraces -1;
    % tracePos --> cellArray
    traceStruct(target_ind).tracePos(target_trace_num) = [];
    % activePause, traceDir --> vector
    traceStruct(target_ind).activePause(target_trace_num)=[];
    traceStruct(target_ind).traceDir(target_trace_num)=[];
    
    
    % These are from measureing RDR
    try
        traceStruct(target_ind).RDR(target_trace_num,:) = [];
        traceStruct(target_ind).RDRFit(target_trace_num,:) = [];
        % uvRange --> cellArray (matrix of cells containing vectors)
        traceStruct(target_ind).uvRange(target_trace_num,:) = [];
    end
    
    % These are from identifying activeSS
    try
        traceStruct(target_ind).TC(target_trace_num)=[];
        traceStruct(target_ind).TC1(target_trace_num)=[];
    end
    try
        % activated
        vec = traceStruct(target_ind).ssActivated;
        if any(vec == target_trace_num)
            vec(vec == target_trace_num) = 0;
            vec(vec > target_trace_num) =  vec(vec > target_trace_num)-1;
            traceStruct(target_ind).ssActivated = vec;
        end
    end
    
    
    % They do exist sometime in history   
    try
        traceStruct(target_ind).TQ1(target_trace_num)=[];
        traceStruct(target_ind).TQ(target_trace_num)=[];
    end
    try
        % traceRDR --> matrix (rows of vectors)
        traceStruct(target_ind).traceRDR(target_trace_num,:) = [];
        traceStruct(target_ind).traceRDRFit(target_trace_num,:) = [];
        traceStruct(target_ind).traceStrain(target_trace_num,:) = [];
        
        traceStruct(target_ind).RDA(target_trace_num,:) = [];
        traceStruct(target_ind).traceDirDeformed(target_trace_num,:) = [];
        traceStruct(target_ind).RDR_adjusted(target_trace_num,:) = [];
        traceStruct(target_ind).RDR_adjustedFit(target_trace_num,:) = [];
        traceStruct(target_ind).RDA_adjusted(target_trace_num,:) = [];
        
    end

    
end

end
