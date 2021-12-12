% These are codes to check if the manually labeled traces are correct.
%% (1) I find that some of the manual labeled twin systems are not the identified active twin systems. 
iE = 5;
v = [];
for iS = 1:length(struCell{iE})
    ID_current = struCell{iE}(iS).gID;
    ind = find(gID==ID_current);
    
    euler = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
    if (1==eulerAligned)
        % g = euler_to_transformation(euler,[0,0,0],[0,0,0]);
        [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [0,0,0], [0,0,0], stressTensor, sampleMaterial, 'twin');
    else
        % g = euler_to_transformation(euler,[-90,180,0],[0,0,0]); % setting-2
        [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler, [-90,180,0], [0,0,0], stressTensor, sampleMaterial, 'twin'); % setting-2
    end
    [ssa, c_a, nss, ntwin, ssGroup] = define_SS(sampleMaterial,'twin');
    traceDir = abs_schmid_factor(nss+1:nss+ntwin,3);
            
    trueTS = sum(struCell{iE}(iS).cTrueTwin,1)>0;
    labeledTS = zeros(1,6);
    for iTwin = 1:6
        if ~isempty(struCell{iE}(iS).tGb{iTwin})
            labeledTS(iTwin) = 1;
        end
    end
    v = [v; trueTS(:),labeledTS(:), traceDir, struCell{iE}(iS).gID * ones(6,1), [1;2;3;4;5;6], struCell{iE}(iS).tSF(:)];
%     struCell{iE}(iS).tGb
%     struCell{iE}(iS).tGbPts
%     struCell{iE}(iS).tGbPtsiE
%     struCell{iE}(iS).tGbNormal    
%     
end
ind = (v(:,1)~=v(:,2));
vv = v(ind,:);
open vv

%% [done] After relabel by importing just the drawn lines and redo the 'manual' label using a modified method, check if some fields were the same as old
% 
% tGb = [];
% tGbPts = [];
% tGbPtsiE = [];
% tGbPtNormal = [];
% 
% ttGb = [];
% ttGbPts = [];
% ttGbPtsiE = [];
% ttGbPtNormal = [];
% 
% for iE = 2:5
%     for iS = 1:length(struCell{iE})
%         
% 
%         
%         for itwin = 1:6
%             for igb = 1:length(struCell{iE}(iS).tGb{itwin})
%                 gbl = struCell{iE}(iS).tGb{itwin}(igb);
%                 tGb = [tGb; gbl ];
%                
%                 for ipt = 1:length(struCell{iE}(iS).tGbPtsiE{itwin}{igb})
%                     tGbPts = [tGbPts; struCell{iE}(iS).tGbPts{itwin}{igb}(ipt,:)];
%                     tGbPtsiE = [tGbPtsiE; struCell{iE}(iS).tGbPtsiE{itwin}{igb}(ipt)];
%                     tGbPtNormal = [tGbPtNormal; struCell{iE}(iS).tGbNormal{itwin}{igb}(ipt)];
%                     
%                 end
%             end
%         end
%         
%         for itwin = 1:6
%             for igb = 1:length(struRef{iE}(iS).tGb{itwin})
%                 gbl = struRef{iE}(iS).tGb{itwin}(igb);
%                 ttGb = [ttGb; gbl ];
%                 for ipt = 1:length(struRef{iE}(iS).tGbPtsiE{itwin}{igb})
%                     ttGbPts = [ttGbPts; struRef{iE}(iS).tGbPts{itwin}{igb}(ipt,:)];
%                     ttGbPtsiE = [ttGbPtsiE; struRef{iE}(iS).tGbPtsiE{itwin}{igb}(ipt)];
%                     ttGbPtNormal = [ttGbPtNormal; struRef{iE}(iS).tGbNormal{itwin}{igb}(ipt)];
%                     
%                 end
%             end
%         end
%         
%         if length(tGbPts)~=length(ttGbPts)
%             [iE,iS]
%             error('rr');
%         end
%         
%     end
% end


%%








