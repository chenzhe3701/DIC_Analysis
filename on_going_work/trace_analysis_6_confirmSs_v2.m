% Zhe Chen, 2016-1-18
% Confirm identified SS, and label slip transfer pairs.
%
% Look at for each Grain
% (1) plot all Traces: label identified SS,
% (2) plot measured RDR with strain.
% (3) Table all SS, SF, and their theoretical RDR
%
% Display existing transferPairs.
% Input ??? if to add transferPair:
% (1)neighbor grain ID #G,
% then update plot with neighbor's labeled traces and trace #.
% then input (2) trace #T in this grain, (3) trace #nT in neighbor
%
% Calculate mPrime between traces
% new Table/quantity mPrime between possible SS in grain 1 and grain 2
% (check direction for possible SS's)
%
% If necessary ???, change identified SS #S of trace #T:
% (Input: (1)Trace #T, (2)new SS #newS, '0' if no ss match)
% So, old SS #S = find(traceStruct().ssActivated == T)
% ssActivated(S)<-0, ssActivated(newS)<-T.
% TC(T)<-0, TC1(T)<-0.
% TQ <- (-1), TQ1 <- 0
% label 'TQ' with -1, so edit these immediately later.
%
% chenzhe, 2017-06-08.  Review code.  Update function.
% Currently, do not consider the field 'TQ', 'TQ1'
%
% chenzhe, 2017-07-25, use try(), so no need to run RDR measurement before
% this step.
%
% required functions:
% ---------- level 1 ----------
% trace_analysis_TiMgAl()

addChenFunction;

[traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
traceStruct = load([traceStructPath,'\',traceStructFile]);
traceStruct = traceStruct.traceStruct;

% trace figures' folder
pathFig = [uigetdir('','select dir of the trace figures'),'\'];
figNameList = dir([pathFig,'\*.fig']);
figNameList = struct2cell(figNameList)';
figNameList = figNameList(:,1);

previousData = 'T5_#7_EbsdToSemForTraceAnalysis';   % The data, that you saved, that contains ebsdToSEM data.
load(previousData);
sampleName=sampleName; sampleMaterial=sampleMaterial; stressTensor=stressTensor; % These are some of the prevously loaded settings.

nameFig = '*.fig';

%%
continueTF = true;      % if continue to open new figures

while continueTF
    nameFig  = uigetfile([pathFig,'\',nameFig], 'select trace figure');
    % If selected a figure, continue analyze.
    if nameFig
        % (1) open figure (OR, redraw?) ????????????????????????
        handleFig = openfig([pathFig,'\',nameFig]);
        
        th3 = findobj(handleFig,'tag','t3');    % table tagged 't3' stores grain id
        idThis = th3.Data;
        indThis = find(cell2mat(arrayfun(@(x) (x.ID==idThis), traceStruct, 'uniformoutput',false)));
        ax1 = findobj(handleFig,'tag','ax1');   % axes tagged 'ax1' has strain map
        
        if ~isempty(indThis)&&(traceStruct(indThis).nTraces>0)    % sometimes, the input id cannot be found
            
            % (1) and annotate labeled traces
            % Also, label the identified SS, if have done so
            for iTrace = 1:traceStruct(indThis).nTraces
                x1 = traceStruct(indThis).tracePos{iTrace}(1);
                y1 = traceStruct(indThis).tracePos{iTrace}(3);
                x2 = traceStruct(indThis).tracePos{iTrace}(2);
                y2 = traceStruct(indThis).tracePos{iTrace}(4);
                plot3(ax1,[x1,x2],[y1,y2],[20,20],'linewidth',2,'color','r');
                identifiedSS = find(traceStruct(indThis).ssActivated == iTrace);
                text((x1+x2)/2,(y1+y2)/2,50,['Tr',num2str(iTrace),', ss',num2str(identifiedSS)],'parent',ax1, 'fontsize',14,'color','m');
            end
            
            
            % (2) plot measured RDR
            try
                colors = {'k','r','b','g','y','c'};
                f2 = figure('position',[50,730,600,250]);
                hold on;
                for iTrace = 1:traceStruct(indThis).nTraces
                    try
                        traceStart = traceStruct(indThis).TQ1(iTrace);  % stop at which trace start to be usable
                    catch
                        traceStart = 1;
                    end
                    if traceStart == 0  % if haven't labeled, use stop 1 as start
                        traceStart = 1;
                    end
                    d.x = traceStart:size(traceStruct(indThis).RDR,2);      % 'RDR' is previously 'traceRDR'
                    % for some added trace, RDR has not been measured.
                    try
                        d.y = traceStruct(indThis).RDR(iTrace, traceStart:end);
                    catch
                        d.y = zeros(1,size(d.x,2));
                    end
                    d.y2 = traceStruct(indThis).ssRDR(find(traceStruct(indThis).ssActivated==iTrace));    % RDR for ss
                    if isempty(d.y2)    % if ss not identified for this trace, then cannot find d.y2, so assign '0'
                        d.y2=0;
                    end
                    plot(d.x, d.y, '-o','markersize',8,'color',colors{iTrace});
                    plot(d.x, repmat(d.y2,1,length(d.x)), '--.','color',colors{iTrace},'markersize',22-4*iTrace);
                end
                warning('off','MATLAB:legend:IgnoringExtraEntries');
                legend('Tr1','predicted-1','Tr2','predicted-2','Tr3','predicted-3','Tr4','predicted-4','Tr5','predicted-5','Tr6','predicted-6');
                xlabel('pause number');ylabel('traceRDR');
                hold off;
                warning('on','MATLAB:legend:IgnoringExtraEntries');
                
                
                % (3) AND plot theoretical RDR for all SS
                ind_current = find(idThis == gID);    % an index of row. '_current' rather than this because related to 'gID3D'
                phi1_current = gPhi1(ind_current);
                phi_current = gPhi(ind_current);
                phi2_current = gPhi2(ind_current);
                [sf_mat, ~, burgersXY] = trace_analysis_TiMgAl([phi1_current,phi_current,phi2_current],[-90,180,0],[180,0,0],stressTensor,sampleMaterial,'noTwin'); %--------------- your sample is flipped !!!!
                f3 = figure('position',[1250,200,350,550]);
                uitable(f3, 'position',[20,20,300,480],...
                    'data', [sf_mat(:,1:2),burgersXY(:,3)],...
                    'ColumnName',{'ss#','SF','RDR'});
            end
            % no transfer pair labeled.   ---------------------------------------------------------------------------------------------
            %             % (4) Operations: input transfer pair? correct identified SS?
            %             if ~isfield(traceStruct(indThis),'transferPair')
            %                 traceStruct(indThis).transferPair=[];
            %             end
            %             disp('-------------------------------------------------------------');   % display existing transferPair
            %             disp(['for grain: grain-ID = ',num2str(idThis),' ,existing transferPair:']);
            %             disp(['grain_1, grain_2, trace_1, trace_2, stremgth']);
            %             disp(traceStruct(indThis).transferPair);
            %
            %
            %
            %
            %             % (4.1) plot existing transfer pair
            %             for ii=1:size(traceStruct(indThis).transferPair,1)
            %                 idNeighbor = traceStruct(indThis).transferPair(ii,2);
            %                 indNeighbor = find(cell2mat(arrayfun(@(x) (x.ID==idNeighbor), traceStruct, 'uniformoutput',false)));
            %                 if ~isempty(indNeighbor)    % sometimes, the input neighbor cannot be found
            %                     iTrace = traceStruct(indThis).transferPair(ii,4);
            %                     x1 = traceStruct(indNeighbor).tracePos{iTrace}(1);
            %                     y1 = traceStruct(indNeighbor).tracePos{iTrace}(3);
            %                     x2 = traceStruct(indNeighbor).tracePos{iTrace}(2);
            %                     y2 = traceStruct(indNeighbor).tracePos{iTrace}(4);
            %                     plot3(ax1,[x1,x2],[y1,y2],[20,20],'linewidth',2,'color','r');
            %                     identifiedSS = find(traceStruct(indNeighbor).ssActivated == iTrace);
            %                     text((x1+x2)/2,(y1+y2)/2,50,['Tr',num2str(iTrace),', ss',num2str(identifiedSS)],'parent',ax1, 'fontsize',14,'color',[.95, .95, .95]);
            %                 end
            %             end
            %
            %
            %
            %             % (4.2) choose if need to input transfer pair
            %             inputTransferPair = 1;  % initiate, inputOperation = YES
            %             while inputTransferPair == 1
            %
            %
            %
            %                 % input if need to input transfer pair
            %                 prompt = ['Input transfer pair: Y/N ?', char(10)];
            %                 str = input(prompt, 's');
            %
            %                 if strcmpi(str, 'Y')
            %                     % input neighbor ID
            %                     prompt = ['Input neighbor grain ID', char(10)];
            %                     str1 = input(prompt, 's');
            %                     str1 = str2num(str1);
            %
            %                     % (4.2.1) update plot with neighbor grains trace
            %                     idNeighbor = str1;
            %                     indNeighbor = find(cell2mat(arrayfun(@(x) (x.ID==idNeighbor), traceStruct, 'uniformoutput',false)));
            %                     if (~isempty(indNeighbor))&&(traceStruct(indNeighbor).nTraces>0)    % sometimes, the input neighbor cannot be found, or neighbor has no trace
            %                         for iTrace = 1:traceStruct(indNeighbor).nTraces
            %                             x1 = traceStruct(indNeighbor).tracePos{iTrace}(1);
            %                             y1 = traceStruct(indNeighbor).tracePos{iTrace}(3);
            %                             x2 = traceStruct(indNeighbor).tracePos{iTrace}(2);
            %                             y2 = traceStruct(indNeighbor).tracePos{iTrace}(4);
            %                             plot3(ax1,[x1,x2],[y1,y2],[20,20],'linewidth',2,'color','r');
            %                             identifiedSS = find(traceStruct(indNeighbor).ssActivated == iTrace);
            %                             text((x1+x2)/2,(y1+y2)/2,50,['Tr',num2str(iTrace),', ss',num2str(identifiedSS)],'parent',ax1, 'fontsize',14,'color',[.95, .95, .95]);
            %                         end
            %
            %                         % (4.2.2) input trace ID in this grain
            %                         prompt = ['Input trace # in this grain', char(10)];
            %                         str2 = input(prompt, 's');
            %                         str2 = str2num(str2);
            %                         % (4.2.3) input trace ID in neighbor grain
            %                         prompt = ['Input trace # in neighbor grain', char(10)];
            %                         str3 = input(prompt, 's');
            %                         str3 = str2num(str3);
            %                         % (4.2.3) input transfer strength
            %                         prompt = ['Input transfer strength, 1 for strong', char(10)];
            %                         str4 = input(prompt, 's');
            %                         if ~strcmp(str4,'1')
            %                             str4 = '0';         % default is 0, weak
            %                         end
            %                         str4 = str2num(str4);
            %
            %
            %                         % (4.2.4) create field 'transferPair'{[this ID, neighborID, this trace, neighbor trace, strength]}
            %                         % also, check if redundant
            %                         tPairToAdd = [idThis, idNeighbor, str2, str3, str4];
            %                         addTransferPairTF = 1;
            %                         for ii = 1:size(traceStruct(indThis).transferPair,1)
            %                             if isequal(tPairToAdd(1:4), traceStruct(indThis).transferPair(ii,1:4));
            %                                 addTransferPairTF = 0;
            %                                 traceStruct(indThis).transferPair(ii,5) = tPairToAdd(5);        % if input pair is the same, do not add, but keep the input strength label.
            %                             end
            %                         end
            %                         if addTransferPairTF
            %                             traceStruct(indThis).transferPair = [traceStruct(indThis).transferPair; tPairToAdd];
            %                         end
            %                         % sort according to this grain ID
            %                         traceStruct(indThis).transferPair = sortrows(traceStruct(indThis).transferPair,2);
            %
            %                         % Also need to add this pair to the neighboring grain !
            %                         tPairToAdd2 = [idNeighbor, idThis, str3, str2, str4];
            %                         addTransferPairTF = 1;
            %                         for ii = 1:size(traceStruct(indNeighbor).transferPair,1)
            %                             if isequal(tPairToAdd2(1:4), traceStruct(indNeighbor).transferPair(ii,1:4));
            %                                 addTransferPairTF = 0;
            %                                 traceStruct(indNeighbor).transferPair(ii,5) = tPairToAdd2(5);
            %                             end
            %                         end
            %                         if addTransferPairTF
            %                             traceStruct(indNeighbor).transferPair = [traceStruct(indNeighbor).transferPair; tPairToAdd2];
            %                         end
            %                         % sort according to this grain ID
            %                         traceStruct(indNeighbor).transferPair = sortrows(traceStruct(indNeighbor).transferPair,2);
            %                     end
            %                 else
            %                     inputTransferPair = 0;
            %                 end
            %
            %
            %             end     % end while
            %
            %
            %
            %
            %             % (4.3) Calculate and display (as Table) mPrime between identified traces
            %             for ii = 1:size(traceStruct(indThis).transferPair,1)
            %                 idNeighbor = traceStruct(indThis).transferPair(ii,2);
            %                 indNeighbor = find(cell2mat(arrayfun(@(x) (x.ID==idNeighbor), traceStruct, 'uniformoutput',false)));
            %                 ind_euler_1 = find(gID3D == idThis);
            %                 ind_euler_2 = find(gID3D == idNeighbor);
            %                 % find (8) closest SS to indexed traces.  If no one is within (10) degrees, choose the closest 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 trace_1 = traceStruct(indThis).transferPair(ii,3);
            %                 trace_1 = find(abs(traceStruct(indThis).ssTraceDir - traceStruct(indThis).traceDir(trace_1)) < 10);
            %                 if isempty(trace_1)
            %                     trace_1 = traceStruct(indThis).transferPair(ii,3);
            %                     [~,trace_1] = sort(abs(traceStruct(indThis).ssTraceDir - traceStruct(indThis).traceDir(trace_1)));
            %                     trace_1 = trace_1(1:8);
            %                 end
            %                 trace_1 = reshape(trace_1,[],1);
            %
            %                 trace_2 = traceStruct(indThis).transferPair(ii,4);
            %                 trace_2 = find(abs(traceStruct(indNeighbor).ssTraceDir - traceStruct(indNeighbor).traceDir(trace_2)) < 10);
            %                 if isempty(trace_2)
            %                     trace_2 = traceStruct(indThis).transferPair(ii,4);
            %                     [~,trace_2] = sort(abs(traceStruct(indNeighbor).ssTraceDir - traceStruct(indNeighbor).traceDir(trace_2)));
            %                     trace_2 = trace_2(1:8);
            %                 end
            %                 trace_2 = reshape(trace_2,1,[]);
            %
            %                 [mPrime,sfG1,sfG2] = calculate_G1G2_mPrime(gPhi13D(ind_euler_1), gPhi3D(ind_euler_1),gPhi23D(ind_euler_1),...
            %                     gPhi13D(ind_euler_2), gPhi3D(ind_euler_2),gPhi23D(ind_euler_2),...
            %                     -90,180,0,0,[1,0,0;0,0,0;0,0,0]);
            %                 mPrime = mPrime(trace_1,trace_2);
            %                 mPrime = [zeros(size(mPrime,1),1),mPrime];
            %                 mPrime = [zeros(1,size(mPrime,2));mPrime];
            %                 mPrime(2:end,1) = sfG1(trace_1);
            %                 mPrime(1,2:end) = sfG2(trace_2);
            %                 % The displayed mPrime has the corresponding Schmid factors in
            %                 % the 1st row and 1st column
            %                 identifiedSS1 = find(traceStruct(indThis).ssActivated == traceStruct(indThis).transferPair(ii,3));
            %                 identifiedSS2 = find(traceStruct(indNeighbor).ssActivated == traceStruct(indThis).transferPair(ii,4));
            %                 if traceStruct(indThis).transferPair(ii,5) == 1
            %                     fn = figure('position',[250,200,800,400],...
            %                         'name',['Strong transferPair, grain: ',num2str(idThis),' ss # ',num2str(identifiedSS1), ' and grain ',num2str(idNeighbor), 'ss # ',num2str(identifiedSS2)] );
            %                 else
            %                     fn = figure('position',[250,200,800,400],...
            %                         'name',['transferPair, grain: ',num2str(idThis),' ss # ',num2str(identifiedSS1), ' and grain ',num2str(idNeighbor), 'ss # ',num2str(identifiedSS2)] );
            %                 end
            %                 uitable(fn, 'position',[20,20,760,355],...
            %                     'data', mPrime,...
            %                     'RowName',['sf,g2->';num2cell(trace_1)],...
            %                     'ColumnName',['sf,g1 v',num2cell(trace_2)]);
            %             end
            
            
            
            
            % (4.4) Input if edit identified slipSystem ?
            pause(1)    % i don't know, but matlab becomes slow.
            correctSS = 1;  % initiate, inputOperation = YES
            while correctSS == 1
                prompt = ['Correct identified active slipSystem: Y/N ?', char(10)];
                str = input(prompt, 's');
                
                if strcmpi(str, 'Y')
                    % input neighbor ID
                    prompt = ['Input trace number', char(10)];
                    str1 = input(prompt, 's');
                    traceNum = str2num(str1);
                    
                    prompt = ['Input new SS.  Enter 0 to make the slip system un-identified', char(10)];
                    str2 = input(prompt, 's');
                    newSsNum = str2num(str2);
                    
                    oldSsNum = find(traceStruct(indThis).ssActivated == traceNum);
                    traceStruct(indThis).ssActivated(oldSsNum) = 0;
                    if newSsNum~=0
                        traceStruct(indThis).ssActivated(newSsNum) = traceNum;
                    end
                    traceStruct(indThis).TC(traceNum) = 0;
                    traceStruct(indThis).TC1(traceNum) = 0;
                    %                     traceStruct(indThis).TQ1(traceNum) = 0;
                    %                     traceStruct(indThis).TQ(traceNum) = -1;
                else
                    correctSS = 0;
                end
            end
            
            try
                close(handleFig);
                close(f2);
                close(f3);
            end
            
        end
        
        close all;
        % find next grain figure name
        try
            nameFig = figNameList{find(strcmpi(nameFig,figNameList))+1};
        catch
            nameFig = figNameList{1};
        end
        
        % Otherwise, stop analyzing
    else
        continueTF = false;
    end
    
end

try
    save('traceStruct_confirmSS_rename_it','traceStruct','-append');
catch
    save('traceStruct_confirmSS_rename_it','traceStruct');
end