% Function add_trace
% run this code in script mode.
% select dir of trace figures, input the pause #, load the traceStructure
% Then draw lines on grains to add/record observed trace.
% 
% The purpose is to record what has been observed, and record
% locations/directions/time of observation.
% The identification of the corresponding/activated slip system will be done in
% other codes.
%
% Zhe Chen, modify on 2015-11-30.
%
% chenzhe, 2017-06-07 review code.  Combine option to remove a trace.
%
% note, required functions in chenFunctions
% ---------- level 1----------
% trace_analysis_3a_removeTrace_v2()


pathFig = [uigetdir('','select dir of the trace figures'),'\'];
figNameList = dir([pathFig,'\*.fig']);
figNameList = struct2cell(figNameList)';
figNameList = figNameList(:,1);
pauseNum = inputdlg('what is the pause # ?');
pauseNum = str2double(pauseNum{1});

[traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
traceStruct = load([traceStructPath,'\',traceStructFile]);
traceStruct = traceStruct.traceStruct;

%%
continueTF = true;      % if continue to open new figures

nameFig = '*.fig';
while continueTF
    nameFig  = uigetfile([pathFig,nameFig], 'select trace figure');
    if nameFig
        handleFig = openfig([pathFig, nameFig]);    % open figure
        
        th3 = findobj(handleFig,'tag','t3');    % table tagged 't3' stores grain id
        idAdd = th3.Data;
        indAdd = find(cell2mat(arrayfun(@(x) (x.ID==idAdd), traceStruct, 'uniformoutput',false)));
        ax1 = findobj(handleFig,'tag','ax1');   % axes tagged 'ax1' has strain map
        
        % if already have trace, or already confirmed, annotate
        for iExistTrace = 1:traceStruct(indAdd).nTraces
            x1 = traceStruct(indAdd).tracePos{iExistTrace}(1);
            y1 = traceStruct(indAdd).tracePos{iExistTrace}(3);
            x2 = traceStruct(indAdd).tracePos{iExistTrace}(2);
            y2 = traceStruct(indAdd).tracePos{iExistTrace}(4);
            plot3(ax1,[x1,x2],[y1,y2],[20,20],'linewidth',2,'color','r');
            identifiedSS = find(traceStruct(indAdd).ssActivated == iExistTrace);
            text((x1+x2)/2,(y1+y2)/2,50,['Tr',num2str(iExistTrace),', ss',num2str(identifiedSS)],'parent',ax1, 'fontsize',14,'color','m');
        end
        
        needToAdd = questdlg('add trace?','select answer','Yes','No','Remove','No');
        switch needToAdd
            case 'Yes'
                % draw trace, get position, add to structure
                nextTrace = true;
                while nextTrace
                    handleImline = imline(ax1);
                    posAdd = wait(handleImline);
                    answer = questdlg('add,redraw,or next grain?','select operation',...
                        'NextTrace','Redraw','AcceptAndNextGrain',...
                        'AcceptAndNextGrain');
                    switch answer
                        case 'NextTrace'
                            directionAdd = atand((posAdd(4)-posAdd(3))/(posAdd(2)-posAdd(1)));
                            if ~isnan(directionAdd)
                                traceNum = traceStruct(indAdd).nTraces + 1;
                                traceStruct(indAdd).nTraces = traceNum;
                                traceStruct(indAdd).tracePos{traceNum} = posAdd;
                                traceStruct(indAdd).traceDir(traceNum) = directionAdd;
                                confirmedPause = inputdlg('confirm active pause?', 'active pause = ', 1, {num2str(pauseNum)});
                                traceStruct(indAdd).activePause(traceNum) = str2num(confirmedPause{1});
                            end
                        case 'Redraw'
                            delete(handleImline);
                        case 'AcceptAndNextGrain'
                            directionAdd = atand((posAdd(4)-posAdd(3))/(posAdd(2)-posAdd(1)));
                            if ~isnan(directionAdd)
                                traceNum = traceStruct(indAdd).nTraces + 1;
                                traceStruct(indAdd).nTraces = traceNum;
                                traceStruct(indAdd).tracePos{traceNum} = posAdd;
                                traceStruct(indAdd).traceDir(traceNum) = directionAdd;
                                confirmedPause = inputdlg('confirm active pause?', 'active pause = ', 1, {num2str(pauseNum)});
                                traceStruct(indAdd).activePause(traceNum) = str2num(confirmedPause{1});
                            end
                            close(handleFig);
                            nextTrace = false;
                            
                            % find next grain figure name
                            try
                                nameFig = figNameList{find(strcmpi(nameFig,figNameList))+1};
                            catch
                                nameFig = figNameList{1};
                            end
                            
                    end
                end
            case 'No'
                close(handleFig);
                % find next grain figure name
                try
                    nameFig = figNameList{find(strcmpi(nameFig,figNameList))+1};
                catch
                    nameFig = figNameList{1};
                end
                
            case 'Remove'
                target_trace_num = inputdlg('which trace number to remove','Input trace number',1);
                target_trace_num = str2double(target_trace_num{1});
                traceStruct = trace_analysis_3a_removeTrace_v2(traceStruct, idAdd, target_trace_num);
                close(handleFig);
        end
    else
        continueTF = false;
    end
    
    % computer crushes often, so save often
    try
        save('traceStruct_for_recover','traceStruct','-append');
    catch
        save('traceStruct_for_recover','traceStruct');
    end
end
%%
answer = questdlg('save structure?');
if strcmpi(answer, 'Yes')
    try
        save('traceStruct_rename_it','traceStruct','-append');
    catch
        save('traceStruct_rename_it','traceStruct');
    end
end

%%
