prompt = {'ID', 'iE:', 'iC:', 'activeSS = [x x x x x x]'};
title_str = 'input ID, iE,iC,activeSS';
dims = [1 50; 1 50; 1 50; 1 50];
defaultInput = {num2str(ids(1)),'0','0','[0 0 0 0 0 0]'};
answer = inputdlg(prompt, title_str, dims, defaultInput);

% If not exist, make one
try
    tNote(1);
catch
    tNote = zeros(1,9);
end

if ~isempty(answer)
    ID_input = str2num(answer{1});
    iE_input = str2num(answer{2});
    iC_input = str2num(answer{3});
    activeSS_input = str2num(answer{4});
    v_input = [ID_input, iE_input, iC_input, activeSS_input];
    
    if (iE_input > 0) && (ID_input > 0) && (iC_input > 0)
        % if modified before, overwrite it.
        for ii=1:size(tNote,1)
            idx = ismember(tNote(:,[1:3]), v_input(:,[1:3]), 'rows');
            if sum(idx)>0
                tNote(idx,:) = zeros(1,9);
            end
        end
        
        tNote = [tNote; v_input];
        
    end
    tNote = unique(tNote,'rows');
    
    saveDataPath = evalin('base','saveDataPath');
    save(fullfile(saveDataPath,'tNote_temp.mat'),'tNote');
end