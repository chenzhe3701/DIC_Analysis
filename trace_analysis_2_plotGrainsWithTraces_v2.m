% ChenZhe, 2015-11-7
% use strain data, EBSD data to plot possible traces for each grain.
% You need to select the (1) the proper folder to save plots
% (2) modify the prefix of the data files.
%
% The purpose is to plot grain, grain #, and strain field, so that you can
% draw lines on these figures to label the observed traces, record their
% directions/time of appearance.
% chenzhe, 2017-06-07.  Review code. Update functions.  Review notes.
% pay attention to FCC.
%
% note, required functions in chenFunctions
% ---------- level 1----------
% grow_boundary()
% trace_analysis_TiMgAl()
% ----> euler_to_transformation()
% ----> crystal_to_cart_ss()
% define_SS()
%
% Fields in structure: 'ID', 'sf', 'ssActivated', 'ssTraceDir', 'ssRDR',
% 'dia', 'e', 'eStd', 'nTraces', 'tracePos', 'activePause'

clear;
addChenFunction;
directory_to_save_plots = uigetdir('','pick a directory to save plots');
directory_DIC = uigetdir('','pick DIC directory, which contains the stitched DIC data for each stop'); 

previousData = 'WE43_T6_C1_EbsdToSemForTraceAnalysis';   % The data, that you saved, that contains ebsdToSEM data. -----------------------------------------
load(previousData);
sampleName=sampleName; sampleMaterial=sampleMaterial; stressTensor=stressTensor; % These are some of the prevously loaded settings.

gIDwithTrace = gID(~isnan(gExx));
% gIDwithTrace = [89,129,135,191,201,210,327,401,422,553];        % WE43 T5 #7, new
gIDwithTrace = [302,151,186,191,296,431,572,1211];            % WE43_T6_C1

% STOP = {'001','002','003','004','005','006','007'};
STOP = {'0','1','2','3','4','5','6','7'};
iE_start = 2;   % elongation levels to analyze. 1-based. ------------------------------------ modify this -----------------
iE_stop = 8;
resReduceRatio = 3;       % to save space, reduce map resolution  %%%%%%%%%%%%%%%%%%% ------------ could modify ----------

% file name prefixes ---------------------------------------------------------------- could modify  --------------- 
f1 = 'T5_#7_stop_';
f2 = '_'; 

notes = struct('atEdge',[],'likeTwin',[],'tooSmall',[]);
% This is a special note for data processing, for specific samples
% useful fields will be 'tooSmall', 'atEdge', 'likeTwin'
% notes = load('T5#7_traceAnalysisNotes');
% 2016-08-15. create a field called 'noTrace'. noTrace = 1 = too small. noTrace = 2 = at edge .

%%
% for iThick = 1:3*resReduceRatio
%     boundaryTF = grow_boundary(boundaryTF);       % grow boundary if necessary  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -------------------------------------------
% end
%%
for iE = iE_start:iE_stop
    clear 'traceStruct';
    strainFile = [directory_DIC,'\',f2,STOP{iE}]; disp(strainFile);            % change the prefix and name of DIC file %%%%%%%%%%%%%%%%%-----------------------------
    

    load(strainFile,'exx');     % Look at exx, but this can be changed in the future. %%%%%%%%%%%%%%%%%%% ------------------------------------------------------------
    directory = [directory_to_save_plots,'\','Trace_plots_stop_',STOP{iE}]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------- change name
    mkdir(directory);     % make new directory to save plots at each stop #
      
    iS = 0;       % index in structure array
    while (iS<length(gIDwithTrace))    % &&(iS<60) % something for debugging 
        % (1) select grain ID, extract data associated with this grain
        iS = iS+1;
        ID_current = gIDwithTrace(iS);              % id of current grain
        ind_current = find(ID_current == gID);    % an index of row
        phi1_current = gPhi1(ind_current);
        phi_current = gPhi(ind_current);
        phi2_current = gPhi2(ind_current);
        
        ID_neighbor = gNeighbors(ind_current,:);
        ID_neighbor = ID_neighbor(ID_neighbor~=0);
        
        % find index range of a small matrix containing the grain of interest
        ind_pool = ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_pool, 1), 1, 'first');
        indC_max = find(sum(ind_pool, 1), 1, 'last');
        indR_min = find(sum(ind_pool, 2), 1, 'first');
        indR_max = find(sum(ind_pool, 2), 1, 'last');
        
        nRow = indR_max - indR_min + 1;
        nColumn = indC_max - indC_min + 1;
        
        e_current = exx(indR_min:indR_max, indC_min:indC_max);  % strain of this region: grain + neighbor. Look at 'exx' strain, but can be changed later %%%%%%%%%%%%%%%
        boundaryTF_current = boundaryTF(indR_min:indR_max, indC_min:indC_max);
        x_current = X(indR_min:indR_max, indC_min:indC_max);
        y_current = Y(indR_min:indR_max, indC_min:indC_max);
        ID_map_current = ID(indR_min:indR_max, indC_min:indC_max);
        e_grain = e_current;
        e_grain(ID_map_current~=ID_current) = nan;  % 'e_grain' is strain of This grain. 'e_current' is strian of this region.
                
        % reduce size  %%% %%%%%%%%%%%%%%%%
        e_current = e_current(1:resReduceRatio:end,1:resReduceRatio:end);
        boundaryTF_current = boundaryTF_current(1:resReduceRatio:end,1:resReduceRatio:end);
        x_current = x_current(1:resReduceRatio:end,1:resReduceRatio:end);
        y_current = y_current(1:resReduceRatio:end,1:resReduceRatio:end);
        ID_map_current = ID_map_current(1:resReduceRatio:end,1:resReduceRatio:end);
        e_grain = e_grain(1:resReduceRatio:end,1:resReduceRatio:end);
                
        % calculate Schmid factor.
        % abs_schmid_factor = [#, SF, angle_XtoY, trace_x_end, trace_y_end]
        % burgersXY = [burgers_X, burgers_Y, ratio].
        % ---------------------------------------- note: Ti or Mg ?  or Al ? ---------------------------------------------------------
        
        [sf_mat, sf_mat_sorted, burgersXY] = trace_analysis_TiMgAl([phi1_current,phi_current,phi2_current],[-90,180,0],[0,0,0],stressTensor,sampleMaterial,'twin');   %--------- M: your sample is flipped
        
        % fill fields
        traceStruct(iS).ID = ID_current;
        if ismember(ID_current, notes.atEdge)
%             traceStruct(iS).notes = 1;
        elseif ismember(ID_current, notes.tooSmall)
%             traceStruct(iS).notes = 2;
        else
%             traceStruct(iS).notes = 0;
            traceStruct(iS).sf = sf_mat(:,2)';
            traceStruct(iS).ssActivated = zeros(size(traceStruct(iS).sf));  % use '0' to hold space for field 'ssActivated'.
            traceStruct(iS).ssTraceDir = sf_mat(:,3)';  % from x-right to y-down coord (CW, that's how line endpoints are read in SEM-plot)
            traceStruct(iS).ssRDR = burgersXY(:,3)';      % changed name 'ssDispRatio' to 'ssRDR'
            
%             [baSF,baSS] = max(sf_mat(1:3,2));
%             [prSF,prSS] = max(sf_mat(4:6,2));   prSS = prSS + 3;
%             [pySF,pySS] = max(sf_mat(7:12,2));  pySS = pySS + 6;
%             [pyCASF,pyCASS] = max(sf_mat(13:18,2)); pyCASS = pyCASS + 12;       % -----------Ti or Mg ---------------------
%             [twinSF,twinSS] = max(sf_mat(19:24,2)); twinSS = twinSS + 18;
%             traceStruct(iS).bppSF = [baSF, prSF, pySF, pyCASF, twinSF];
%             traceStruct(iS).bppSS = [baSS, prSS, pySS, pyCASS, twinSS];
%             traceStruct(iS).bppActivated = zeros(1,5);      % -----------Ti or Mg ---------------------   -
            
            traceStruct(iS).dia = gDiameter(ind_current);
            traceStruct(iS).e = nanmean(e_grain(:));
            traceStruct(iS).eStd = nanstd(e_grain(:));
            
            traceStruct(iS).nTraces = 0;
            traceStruct(iS).tracePos = [];
            traceStruct(iS).activePause = [];
            
            % (2) -------------- plot ---------------------
            f = figure;
            set(f,'position',[50,100,1200,600]);
            % (2.1) Plot grain strain map on left
            ax1 = axes('parent',f,'unit','normalized','position',[0.05, 0.1, 0.4, 0.8]);
            surf(ax1, x_current, y_current, e_current-10,'edgecolor','none');        % surf() can eliminate tag
            cbar = colorbar;
            c1 = quantile(e_grain(~isnan(e_grain)),0.03)-10 - 100*eps;
            c2 = quantile(e_grain(~isnan(e_grain)),0.97)-10 + 100*eps;
%             c1 = nanmean(e_grain(:)) - 2*nanstd(e_grain(:))-10 - 100*eps;
%             c2 = nanmean(e_grain(:)) + 2*nanstd(e_grain(:))-10 + 100*eps;
            if ~isnan(c1+c2)
                caxis(ax1, [c1,c2]);
            end
            cTickLabels = get(cbar, 'ticklabels');      % makes it easier to draw imline later
            for iLabel = 1:length(cTickLabels)
                cTickLabels{iLabel} = num2str(str2num(cTickLabels{iLabel})+10);
            end
            set(cbar,'ticklabels',cTickLabels);
            
            hold on;
            boundaryTF_current(boundaryTF_current==0)=NaN;
            boundaryTF_curent = boundaryTF_current * max(max(e_current(:)),10000);
            surf(x_current, y_current, boundaryTF_current);
            pos_text_y = y_current(ID_map_current == ID_current);
            pos_text_y = nanmean(pos_text_y);
            pos_text_x = x_current(ID_map_current == ID_current);
            pos_text_x = nanmean(pos_text_x);
            pos_text_z = nanmax(c2+1,1);
            text(pos_text_x,pos_text_y,pos_text_z,num2str(ID_current),'fontsize',16);
            
            for ii = 1:length(ID_neighbor)
                ID_nb = ID_neighbor(ii);
                %         ind_nb = find(ID_nb == gID3D);
                pos_text_y = y_current(ID_map_current == ID_nb);
                pos_text_y = nanmean(pos_text_y);
                pos_text_x = x_current(ID_map_current == ID_nb);
                pos_text_x = nanmean(pos_text_x);
                pos_text_z = nanmax(c2+1,1);
                text(pos_text_x,pos_text_y,pos_text_z,num2str(ID_nb),'fontsize',16,'color','r');
            end
            
            axis equal;
            %     [limit_y, limit_x] = size(e_current);
            %     set(ax1,'Ydir','reverse','xLim',[1,limit_x],'yLim',[1,limit_y],'tag','ax1');
            set(ax1,'Ydir','reverse','tag','ax1');
            view(0,90);
            
            % (2.2) Plot traces on the right
            [~,~,nss,nts,ssGroup] = define_SS(sampleMaterial,'twin'); % --------------------------------------------------------- if have twin, change this.
            ssGroup(length(ssGroup)+1:6)=ssGroup(end);
            ax2 = axes('parent',f,'unit','normalized','position',[0.5, 0.25, 0.25, 0.5],'tag','ax2');
            hold on;
            set(ax2, 'xlim', [-1.2, 1.2], 'ylim', [-1.2, 1.2],'ydir','reverse');
            set(ax2 ,'xticklabel','','yticklabel','');
            axis square;
            
            for k = 1:1:nss+nts      % Plot traces
                if k<=ssGroup(1)
                    plot([0,sf_mat(k,4)], [0,sf_mat(k,5)],'r','linewidth',3);
                    text(sf_mat(k,4)*(2+rem(k,17))/17, sf_mat(k,5)*(2+rem(k,17))/17, num2str(k),'color','r','fontsize',12);
                elseif k<=ssGroup(2)
                    plot([0,sf_mat(k,4)], [0,sf_mat(k,5)],'b','linewidth',2);
                    text(sf_mat(k,4)*(2+rem(k,17))/17, sf_mat(k,5)*(2+rem(k,17))/17, num2str(k),'color','b','fontsize',12);
                elseif k<=ssGroup(3)
                    plot([0,sf_mat(k,4)], [0,sf_mat(k,5)],'color',[0,0.5,0],'linewidth',1.5);
                    text(sf_mat(k,4)*(2+rem(k,17))/17, sf_mat(k,5)*(2+rem(k,17))/17, num2str(k),'color',[0,0.5,0],'fontsize',12);
                elseif k<=ssGroup(4)
                    plot([0,sf_mat(k,4)], [0,sf_mat(k,5)],'k','linewidth',1);
                    text(sf_mat(k,4)*(2+rem(k,17))/17, sf_mat(k,5)*(2+rem(k,17))/17, num2str(k),'color','k','fontsize',12);
                elseif k<=ssGroup(5)
                    plot([0,sf_mat(k,4)], [0,sf_mat(k,5)],'c','linewidth',1);
                    text(sf_mat(k,4)*(2+rem(k,17))/17, sf_mat(k,5)*(2+rem(k,17))/17, num2str(k),'color','c','fontsize',12);
                elseif k<=ssGroup(6)
                    plot([0,sf_mat(k,4)], [0,sf_mat(k,5)],'c','linewidth',1);
                    text(sf_mat(k,4)*(2+rem(k,17))/17, sf_mat(k,5)*(2+rem(k,17))/17, num2str(k),'color','m','fontsize',12);
                end
            end
            hold off;
            
            %( 2.3) plot SFs
            t1 = uitable(f,'units','normalized','position',[.8, .05, .18, .9],'Data',sf_mat_sorted(:,1:2), 'ColumnName',{'ssNum','SF'},'ColumnWidth',{50,50},'tag','t1');
            t2 = uitable(f,'Data',[sf_mat(:,1:2),-sf_mat(:,3)],'ColumnName',{'ssNum', 'SF', 'dir'}, 'visible','off','tag','t2');
            t3 = uitable(f,'Data',ID_current,'ColumnName',{'ID'}, 'visible','off','tag','t3');
            % display('schmid factors:');
            % disp(sf_mat);
            
            savefig(f,[directory,'\','trace_grain_',num2str(ID_current),'_stop_',STOP{iE}]);
            saveas(f,[directory,'\','trace_grain_',num2str(ID_current),'_stop_',STOP{iE}],'tif');
            close(f);
            
        end
    end
    
    %traceM = cell2mat(struct2cell(traceStruct)');
    
    % save the structure,
    save([directory,'\','traceStruct_stop_',STOP{iE}],'traceStruct');
    
end



