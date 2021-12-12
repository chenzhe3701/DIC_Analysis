% estimate the grip direction based on xyTrans and xyStage info.
%
% chenzhe, 2017-05-22
%
% chenzhe, 2017-09-20 added note
% (1) This is truly a historic code for  20170430 test.
% (2) It looks at the hdr file for the four corner fovs (could look at all
% fovs).  Analyze the working distance, and stage position change.
% (3) The plot is a stage position change, indicating the axis of sample
% elongation -- this sample elongated along a direction 2.74 degrees ccw
% wrt the horizontal direction.
% (4) You could also look at WD, but looks like this sample was OK.
% (5) The use of this code is very specific.  So I would not do too much
% optimization or generalization.
%
% chenzhe, 2018-02-04, based on m_11_evaluate_grip_align()
% study the working distance, for FEI images.

pathHdr = uigetdir('D:\WE43_T6_C1_insitu_compression\iFast_auto_hdr_file','choose the hdr folder, all hdr in one folder');
pathDic = uigetdir('D:\WE43_T6_C1_insitu_compression\allFov','choose the dic folder, all dic in one folder');
saveDataPath = [uigetdir('','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];
if ~strcmpi(pathHdr(end),'\')
    pathHdr = [pathHdr,'\'];
end
if ~strcmpi(pathDic(end),'\')
    pathDic = [pathDic,'\'];
end

f1 = 'WE43_T6_C1_s';
ne = 7;     % don't change e,r,c to E,R,C
nr = 10;
nc = 18;
micronPerPixel = 540/6144;
for ie = 0:ne
    for ir = 0:nr
        for ic = 0:nc
            iE = ie + 1;
            iR = ir + 1;
            iC = ic + 1;
            
            fName = [f1,num2str(ie),'_r',num2str(ir),'c',num2str(ic)];
            fNameHdr = [fName,'.hdr'];
            fNameDic = [fName,'.mat'];
            
            WD(iR,iC,iE) = get_hdr_field_value([pathHdr,fNameHdr],'WD');
            StageX(iR,iC,iE) = get_hdr_field_value([pathHdr,fNameHdr],'StageX');
            StageY(iR,iC,iE) = get_hdr_field_value([pathHdr,fNameHdr],'StageY');
            
            WD(iR,iC,iE) = WD(iR,iC,iE) * 1e6;  % convert to micron
            StageX(iR,iC,iE) = StageX(iR,iC,iE) * 1e6; % make it to micron
            StageY(iR,iC,iE) = -StageY(iR,iC,iE) * 1e6;    % make it the same as image coordinate
            
            load([pathDic,fNameDic],'u','v');
            [nR,nC] = size(u);
            tR = round(nR/2);
            tC = round(nC/2);
            ImageU(iR,iC,iE) = u(tR,tC)*micronPerPixel;
            ImageV(iR,iC,iE) = v(tR,tC)*micronPerPixel;
            
            dW(iR,iC,iE) = WD(iR,iC,iE) - WD(iR,iC,1);
            dX(iR,iC,iE) = StageX(iR,iC,iE) - StageX(iR,iC,1) + ImageU(iR,iC,iE) - ImageU(iR,iC,1);
            dY(iR,iC,iE) = StageY(iR,iC,iE) - StageY(iR,iC,1) + ImageV(iR,iC,iE) - ImageV(iR,iC,1);
        end
    end
end
save([saveDataPath,'stageDataInTest.mat'],'WD','StageX','StageY','dW','dX','dY');
%% This is to estimate and plot sample rotation.  However, note it assumes that there is no slip between sample and grip, so that imaging position can represent grip move direction 
% figure; hold on;set(gca,'ydir','reverse');axis equal;
% colors = parula(ne+1);
% data = [];
% for ie = 0:ne
%     for ir = [0:nr]
%         for ic = [0,nc]
%             iE = ie + 1;
%             iR = ir + 1;
%             iC = ic + 1;
%             plot(dX(iR,iC,iE),dY(iR,iC,iE),'o','color',colors(iE,:));
%             data = [data;dX(iR,iC,iE),dY(iR,iC,iE)];
%         end
%     end
% end
% model = fitlm(data(:,1),data(:,2))
% angle = atand(model.Coefficients.Estimate(2))
% yOffSet = micronPerPixel*69000*model.Coefficients.Estimate(2)

%%

figure; hold on;set(gca,'ydir','reverse','zdir','reverse');
% axis equal;
colors = parula(ne+1);
p = [];
tilt = [];
for ie = 0:ne
    iE = ie+1;
    for ir = 1:nr-1
       for ic = 1:nc-1
          xx = StageX(ir:ir+1,ic:ic+1,iE);
          yy = StageY(ir:ir+1,ic:ic+1,iE);
          zz = WD(ir:ir+1,ic:ic+1,iE);     % make negative
          fill3(xx([1,2,4,3]),yy([1,2,4,3]),zz([1,2,4,3]),colors(iE,:),'facealpha',0.8);
       end
    end
    xx = StageX(:,:,iE);
    yy = StageY(:,:,iE);
    zz = WD(:,:,iE);
    model = fit([xx(:),yy(:)],zz(:),'poly11')
    p = [p; model.p00, model.p10, model.p01];
    tilt = [tilt; atand(model.p10), atand(model.p01)];
end
xlabel('stage x, micron');
ylabel('stage y, micron');
zlabel('working distance, micron');



