% ZheChen update, 2016-4-26
% Also measure the range of u and v of the data selected.
% Also measure the rotation of the slip trace, i.e., grain.
%
% ZheChen, 2015-12-6
% For each strain level, measure displacement ratio of the labeled traces.
% Also, measure the strain level of each labeled trace at each interval.
%
% loop-1: strain files
% loop-2: each grain
% loop-3: each trace
%
% chenzhe 2017-07-07.  remove old method.  Review code.
%
% note, required functions in chenFunctions
% ---------- level 1----------
% grids_covered_by_line()
%
% new fields after this function: 'RDR', 'RDRFit', 'uvRange'

addChenFunction;

STOP = {'001','002','003','004','005','006','007'};
iE_start = 3;   % elongation levels to analyze. 1-based. ------------------------------------ modify this -----------------
iE_stop = length(STOP); %% ------------------------- default but could modify this -----------------
% file name prefixes ---------------------------------------------------------------- could modify  ---------------
f1 = '';    %'20170430_ts5Al_02_e';
f2 = '_';

directory_DIC = uigetdir('','pick DIC directory');

[traceStructFile, traceStructPath] = uigetfile('','select traceStructure');
traceStruct = load([traceStructPath,'\',traceStructFile]);
traceStruct = traceStruct.traceStruct;
nS = length(traceStruct);

previousData = 'T5_#7_EbsdToSemForTraceAnalysis';   % The data, that you saved, that contains ebsdToSEM data.
load(previousData);

%%
for iE = iE_start:iE_stop
    
    strainFile = [directory_DIC,'\',f1,STOP{iE}]; disp(strainFile);    % note/change the prefix of DIC file %%%%%%%%%%%%%%%%%%%%
    load(strainFile,'x','y','u','v','exx','sigma');     % This can be changed in the future. %%%%%%%%%%%%%%%%%%%
    inds = (sigma==-1);
    u(inds) = nan;
    v(inds) = nan;
    exx(inds) = nan;
    
    for iS = 1:nS
        disp(['grain: ',num2str(traceStruct(iS).ID)]);
        for iTrace = 1:traceStruct(iS).nTraces
            
            dataRange = 5;      % currently, data range is 11 points in total ----------- this could/should be modified and is related to window/step size. --------------
            
            % method (1) is omitted. -- the horizontal line method
            % method (2) take another line that is perpendicular to the trace line, then take data along this perpendicular line
            
            pt1 = traceStruct(iS).tracePos{iTrace}(1,:);
            pt2 = traceStruct(iS).tracePos{iTrace}(2,:);
            [pt,ind,indR,indC] = grids_covered_by_line(x, y, pt1, pt2);     % pt=[x,y] position, and inds of the data points on the slip trace line.
            
            %             [~,ind1_r] = min(abs(y(:,1)-pt1(2)));
            %             [~,ind1_c] = min(abs(x(1,:)-pt1(1)));
            %             [~,ind2_r] = min(abs(y(:,1)-pt2(2)));
            %             [~,ind2_c] = min(abs(x(1,:)-pt2(1)));
            
%             % Here we need to find the direction of the deformed slip trace. This is not so straight forward, because it could intersecting a stronger slip trace line of another slip system.
%             % Therefore we could take two points, located not so far away on this slip trace line -- so hopefully there are no intersecting lines between them.  Find the direction, and do this many times ...
%             for ii = 1:size(pt,1)-11
%                 tAngle_deformed(ii) = atand( (y(indR(ii+10),indC(ii+10))+v(indR(ii+10),indC(ii+10))-y(indR(ii),indC(ii))-v(indR(ii),indC(ii)))/(x(indR(ii+10),indC(ii+10))+u(indR(ii+10),indC(ii+10))-x(indR(ii),indC(ii))-u(indR(ii),indC(ii))) );
%             end
%             tAngle_deformed = nanmean(tAngle_deformed);
%             %             tAngle_deformed = atand( (pt2(2)+v(ind2_r,ind2_c)-pt1(2)-v(ind1_r,ind1_c))/(pt2(1)+u(ind2_r,ind2_c)-pt1(1)-u(ind1_r,ind1_c)) );
%             %             % sometimes, the u or v is nan.
%             %             if isnan(tAngle_deformed)
%             %                 tAngle_deformed = atand( (pt2(2)+v(ind2_r-2:ind2_r+2,ind2_c-2:ind2_c+2)-pt1(2)-v(ind1_r-2:ind1_r+2,ind1_c-2:ind1_c+2))./(pt2(1)+u(ind2_r-2:ind2_r+2,ind2_c-2:ind2_c+2)-pt1(1)-u(ind1_r-2:ind1_r+2,ind1_c-2:ind1_c+2)) );
%             %                 tAngle_deformed = nanmean(tAngle_deformed(:));
%             %             end
            
            
            uv = [];
            uvAdjusted = [];
            uv_range = [];
            tAngle = atand((pt2(2)-pt1(2))/(pt2(1)-pt1(1)));   % Angle of slip trace line. assume line is not verticle.
            pAngle = tAngle + 90;   % angle of the perpendicular line of slip trace line.
            indRange = ceil(max(abs([sind(pAngle),cosd(pAngle)]))*dataRange);   % size of the window (unit in pixel, or index) of local data points considered
            %             tAngle_change = tAngle_deformed - tAngle;
            
            for ii = 1:length(ind)
                indr = indR(ii);    % index of the current data point on slip trace line
                indc = indC(ii);
                
                xLocal = x(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
                yLocal = y(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
                uLocal = u(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
                vLocal = v(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
                exxLocal = exx(indr-indRange:indr+indRange,indc-indRange:indc+indRange);
                
                [ptLocal,indLocal,indrLocal,indcLocal] = grids_covered_by_line(xLocal, yLocal, pt(ii,:)-[1000,1000*tand(pAngle)], pt(ii,:)+[1000,1000*tand(pAngle)]);
                
                uvLocal = [uLocal(indLocal),vLocal(indLocal)]';
                row_mean = nanmean(uvLocal,2);
                uvLocal = uvLocal - repmat(row_mean,1,size(uvLocal,2));      % displacements u,v, and post-deformed locations x_p,y_p are all centered.
                uv = [uv,uvLocal];
                uv_range = [uv_range,max(uvLocal,[],2)-min(uvLocal,[],2)];
                
                % THe artificial relative displacement is determined by (1) the amount of ration, and (2) the position of these data points from which the RDR is measured
                %                 g = [cosd(-tAngle_change), cosd(-tAngle_change - 90);
                %                     cosd(-tAngle_change + 90), cosd(-tAngle_change)];
                %
                %                 x_deformed_local = xLocal + uLocal;
                %                 y_deformed_local = yLocal + vLocal;
                %                 xyLocal = [xLocal(indLocal),yLocal(indLocal)]';     % reference position
                %                 xy_deformed_local = [x_deformed_local(indLocal),y_deformed_local(indLocal)]';  % deformed position
                %                 xy_deformed_adjusted = g'* xy_deformed_local;       % rotated the measured x_p,y_p position back, eliminating the rotation effect
                %                 uvLocalAdjusted = xy_deformed_adjusted - xyLocal;
                %                 row_mean = nanmean(uvLocalAdjusted,2);
                %                 uvLocalAdjusted = uvLocalAdjusted - repmat(row_mean,1,size(uvLocalAdjusted,2));
                %                 uvAdjusted = [uvAdjusted,uvLocalAdjusted];
                
                %             % this is used in debugging, try to record a trace line locally.  Reuse
                %             % it with modification.
                %             lineLocal = [];
                %             for jj = 95:105
                %                 lineLocal = [lineLocal; x(indR(jj),indC(jj)),y(indR(jj),indC(jj))];
                %             end
                %             % end of debug part
            end
            
            % calculate RDR and RDA
            %             figure;hold on;
            %             scatter(uv(2,:),uv(1,:));
            %             h = lsline;
            %             slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
            ok = ~sum(isnan(uv),1);
            uv = uv(:,ok);
            coef = corrcoef(uv(1,:)',uv(2,:)');
            %             traceStruct(iS).RDRCorr(iTrace,iFile) = coef(2);
            t1 = fitlm(uv(2,:),uv(1,:));
            traceStruct(iS).RDR(iTrace,iE) = t1.Coefficients{2,1};
            %             traceStruct(iS).RDRError(iTrace,iFile) = t1.Coefficients{2,2};
            %             traceStruct(iS).RDRRMSE(iTrace,iFile) = t1.RMSE;
            traceStruct(iS).RDRFit(iTrace,iE) = t1.Rsquared.Ordinary;
            %             traceStruct(iS).RDA(iTrace,iFile) = acotd(t1.Coefficients{2,1});
            %
            %             %             traceStruct(iS).RDR(iTrace,iFile) = slope;
            %             %             traceStruct(iS).RDA(iTrace,iFile) = acotd(slope);
            %             traceStruct(iS).traceDirDeformed(iTrace,iFile) = tAngle_deformed;
            traceStruct(iS).uvRange{iTrace,iE} = quantile(uv_range,0.99,2)-quantile(uv_range,0.01,2);
            %             %             close all;
            %             %             clear h;
            %
            %             % calculate adjusted RDR and RDA
            %             %             figure;hold on;
            %             %             scatter(uvAdjusted(2,:),uvAdjusted(1,:));
            %             %             h = lsline;
            %             %             slope = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
            %             ok = ~sum(isnan(uvAdjusted),1);
            %             uvAdjusted = uvAdjusted(:,ok);
            %             coef = corrcoef(uvAdjusted(1,:)',uvAdjusted(2,:)');
            % %             traceStruct(iS).RDR_adjustedCorr(iTrace,iFile) = coef(2);
            %             t2 = fitlm(uvAdjusted(2,:),uvAdjusted(1,:));
            %             traceStruct(iS).RDR_adjusted(iTrace,iFile) = t2.Coefficients{2,1};
            % %             traceStruct(iS).RDR_adjustedError(iTrace,iFile) = t2.Coefficients{2,2};
            % %             traceStruct(iS).RDR_adjustedRMSE(iTrace,iFile) = t2.RMSE;
            %             traceStruct(iS).RDR_adjustedFit(iTrace,iFile) = t2.Rsquared.Ordinary;
            %             traceStruct(iS).RDA_adjusted(iTrace,iFile) = acotd(t2.Coefficients{2,1});
            %             %             traceStruct(iS).RDR_adjusted(iTrace,iFile) = slope;
            %             %             traceStruct(iS).RDA_adjusted(iTrace,iFile) = acotd(slope);
            %             %             close all;
            %             %             clear h;
        end
    end
    try
        save([traceStructPath,'\traceStruct_RDR_measured_rename_it'],'traceStruct','-append');
        disp('save append');
    catch
        save([traceStructPath,'\traceStruct_RDR_measured_rename_it'],'traceStruct');
        disp('saved');
    end
end

%%
answer = questdlg('save structure?');
if strcmpi(answer, 'Yes')
    traceStruct = traceStruct;
    try
        save([traceStructPath,'\traceStruct_RDR_measured_rename_it'],'traceStruct','-append');
        disp('save append');
    catch
        save([traceStructPath,'\traceStruct_RDR_measured_rename_it'],'traceStruct');
    end
end







