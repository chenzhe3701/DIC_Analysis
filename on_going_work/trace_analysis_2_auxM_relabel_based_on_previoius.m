% This code used previous confirmed ground truth data as a reference to
% fast process data run by updated codes.
% The codes are being cleaned up after 2nd paper acceptance.
%
% chenzhe, 2018-08-27
%
% find groud truth data, use it to correct twins.

ref_file = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\WE43_T6_C1_s5_cluster_to_twin_result.mat';
stru_ref = load(ref_file,'stru');
stru_ref = stru_ref.stru;

%%
ii = 0;
%%
continueTF = true;
try 
    hh.delete;
    hhh.delete;
catch
end
while continueTF && (ii<length(stru))
    ii = ii + 1
    jj = find(arrayfun(@(x) x.gID==stru(ii).gID, stru_ref));
    % make sure operation is on the same grain
    
    
    % if labeld as twin, check if need to disable
    % otherwise, check if need to enable
    if any(stru(ii).c2t)
        continueTF = false;
        disp('check if need to disable');
        disp(['grain: ',num2str(ii)]);
        disp(stru(ii).c2t');
        try
            disp(stru_ref(jj).trueTwin');
        catch
            disp('no corresponding grain in stru_ref');
        end
    elseif ~isempty(jj) && any(stru_ref(jj).trueTwin)
        continueTF = false;
        disp('check if need to enable');
        disp(['grain: ',num2str(ii)]);
        disp(stru(ii).c2t');
        disp(stru_ref(jj).trueTwin')
    end
    
    if continueTF == false
        ind = find(gID == stru(ii).gID);
        hh = imellipse(at,[gCenterX(ind)-500,gCenterY(ind)-500,1000,1000]);
        hhh = imellipse(at,[gCenterX(ind)-2000,gCenterY(ind)-2000,4000,4000]);
    end
    
end