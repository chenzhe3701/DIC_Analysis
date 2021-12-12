% ssa = define_SS(str1, str2); give material name in str1, such as 'Ti', 'pyii'.
% If want to include twin, use str2 = 'twin'.
% Return ssa which is slip systems.
% ssa(1,:,#) = slip plane for slip system #
% ssa(2,:,#) = slip direction for slip system #
%
% added pyii<c+a> slip on Jan-29-2016
%
% Zhe Chen, 2016-1-29
% 
% chenzhe, 2017-06-07. Modify a little bit. Add FCC 'Al'
% add nss, ntwin, ssGroup as output.
% ssGroup is just to aid plotting.
% haven't add twin system for FCC.
%
% chenzhe, 2019-09-11. add contraction twin 'ctwin' to Mg.
% did not check if correct.
% 
% chenzhe, add note 2020-03-25
% choice for str2 = {'twin','ctwin','pyii','whatever=notwin'}

function [ssa, c_a, nss, ntwin, ssGroup] = define_SS(str1,str2)

if any(strcmpi(str1,{'Ti','Titanium','Zr','Zirconium'}))
    c_a = 1.59;
    nss = 24;
    ntwin = 0;
    ssGroup = [3,6,12,24];
    %  basal <a>-glide
    ssa(:,:,1) = [0 0 0 1; 2 -1 -1 0];
    ssa(:,:,2) = [0 0 0 1; -1 2 -1 0];
    ssa(:,:,3) = [0 0 0 1; -1 -1 2 0];
    %  prism <a>-glide
    ssa(:,:,4) = [0 1 -1 0; 2 -1 -1 0];
    ssa(:,:,5) = [1 0 -1 0; -1 2 -1 0];
    ssa(:,:,6) = [-1 1 0 0; -1 -1 2 0];
    %  pyramidal <a>-glide
    ssa(:,:,7) = [0 1 -1 1; 2 -1 -1 0];
    ssa(:,:,8) = [1 0 -1 1; -1 2 -1 0];
    ssa(:,:,9) = [-1 1 0 1; -1 -1 2 0];
    ssa(:,:,10) = [0 -1 1 1; 2 -1 -1 0];
    ssa(:,:,11) = [-1 0 1 1; -1 2 -1 0];
    ssa(:,:,12) = [1 -1 0 1; -1 -1 2 0];
    %  pyramidal <c+a>-glide
    % % This is not consistent with pyramidal.  I changed it on 2016-Jan20
    %     ssa(:,:,13) = [1 0 -1 1; 2 -1 -1 -3];
    %     ssa(:,:,14) = [1 0 -1 1; 1 1 -2 -3];
    %     ssa(:,:,15) = [0 1 -1 1; 1 1 -2 -3];
    %     ssa(:,:,16) = [0 1 -1 1; -1 2 -1 -3];
    %     ssa(:,:,17) = [-1 1 0 1; -1 2 -1 -3];
    %     ssa(:,:,18) = [-1 1 0 1; -2 1 1 -3];
    %     ssa(:,:,19) = [-1 0 1 1; -2 1 1 -3];
    %     ssa(:,:,20) = [-1 0 1 1; -1 -1 2 -3];
    %     ssa(:,:,21) = [0 -1 1 1; -1 -1 2 -3];
    %     ssa(:,:,22) = [0 -1 1 1; 1 -2 1 -3];
    %     ssa(:,:,23) = [1 -1 0 1; 1 -2 1 -3];
    %     ssa(:,:,24) = [1 -1 0 1; 2 -1 -1 -3];
    ssa(:,:,13) = [0 1 -1 1; -1 -1 2 3];
    ssa(:,:,14) = [0 1 -1 1; 1 -2 1 3];
    ssa(:,:,15) = [1 0 -1 1; -2 1 1 3];
    ssa(:,:,16) = [1 0 -1 1; -1 -1 2 3];
    ssa(:,:,17) = [-1 1 0 1; 1 -2 1 3];
    ssa(:,:,18) = [-1 1 0 1; 2 -1 -1 3];
    ssa(:,:,19) = [0 -1 1 1; 1 1 -2 3];
    ssa(:,:,20) = [0 -1 1 1; -1 2 -1 3];
    ssa(:,:,21) = [-1 0 1 1; 2 -1 -1 3];
    ssa(:,:,22) = [-1 0 1 1; 1 1 -2 3];
    ssa(:,:,23) = [1 -1 0 1; -1 2 -1 3];
    ssa(:,:,24) = [1 -1 0 1; -2 1 1 3];
    
    % pyramidal ii <c+a>
    if strcmpi(str2,'pyii')
        nss = 30;
        ssGroup = [ssGroup,30];
        ssa(:,:,25) = [1 1 -2 2; -1 -1 2 3];
        ssa(:,:,26) = [-1 2 -1 2; 1 -2 1 3];
        ssa(:,:,27) = [-2 1 1 2; 2 -1 -1 3];
        ssa(:,:,28) = [-1 -1 2 2; 1 1 -2 3];
        ssa(:,:,29) = [1 -2 1 2; -1 2 -1 3];
        ssa(:,:,30) = [2 -1 -1 2; -2 1 1 3];
    end
    
    % twin system
    % T1: {10-12}<-1011>
    % T2: {11-21><-1-126>
    % C1: {11-22}<11-2-3>
    % C2: {10-11}<10-12>
    if strcmpi(str2,'twin')
        ntwin = 6;
        ssGroup = [ssGroup,ssGroup(end)+6];
        % extension twin type-1
        ssa(:,:,25) = [1 0 -1 2; -1 0 1 1];
        ssa(:,:,26) = [0 1 -1 2; 0 -1 1 1];
        ssa(:,:,27) = [-1 1 0 2; 1 -1 0 1];
        ssa(:,:,28) = [-1 0 1 2; 1 0 -1 1];
        ssa(:,:,29) = [0 -1 1 2; 0 1 -1 1];
        ssa(:,:,30) = [1 -1 0 2; -1 1 0 1];
    end
    if strcmpi(str2,'ctwin')
        % t1 + c1
        ntwin = 12;
        ssGroup = [ssGroup,ssGroup(end)+6,ssGroup(end)+12];
        ssa(:,:,25) = [1 0 -1 2; -1 0 1 1];
        ssa(:,:,26) = [0 1 -1 2; 0 -1 1 1];
        ssa(:,:,27) = [-1 1 0 2; 1 -1 0 1];
        ssa(:,:,28) = [-1 0 1 2; 1 0 -1 1];
        ssa(:,:,29) = [0 -1 1 2; 0 1 -1 1];
        ssa(:,:,30) = [1 -1 0 2; -1 1 0 1];
        
        ssa(:,:,31) = [1 1 -2 2; 1 1 -2 -3];
        ssa(:,:,32) = [-1 -2 -1 2; -1 -2 -1 -3];
        ssa(:,:,33) = [-2 1 1 2; -2 1 1 -3];
        ssa(:,:,34) = [-1 -1 2 2; -1 -1 2 -3];
        ssa(:,:,35) = [1 -2 1 2; 1 -2 1 -3];
        ssa(:,:,36) = [2 -1 -1 2; 2 -1 -1 -3];
    end
    
elseif any(strcmpi(str1,{'Mg','Magnesium'}))
    c_a = 1.62;
    nss = 18;
    ntwin = 0;
    ssGroup = [3,6,12,18];
    %  basal <a>-glide
    ssa(:,:,1) = [0 0 0 1; 2 -1 -1 0];
    ssa(:,:,2) = [0 0 0 1; -1 2 -1 0];
    ssa(:,:,3) = [0 0 0 1; -1 -1 2 0];
    %  prism <a>-glide
    ssa(:,:,4) = [0 1 -1 0; 2 -1 -1 0];
    ssa(:,:,5) = [1 0 -1 0; -1 2 -1 0];
    ssa(:,:,6) = [-1 1 0 0; -1 -1 2 0];
    %  pyramidal <a>-glide
    ssa(:,:,7) = [0 1 -1 1; 2 -1 -1 0];
    ssa(:,:,8) = [1 0 -1 1; -1 2 -1 0];
    ssa(:,:,9) = [-1 1 0 1; -1 -1 2 0];
    ssa(:,:,10) = [0 -1 1 1; 2 -1 -1 0];
    ssa(:,:,11) = [-1 0 1 1; -1 2 -1 0];
    ssa(:,:,12) = [1 -1 0 1; -1 -1 2 0];
    %  pyramidal <c+a>-glide, 2nd order
    ssa(:,:,13) = [1 1 -2 2; -1 -1 2 3];
    ssa(:,:,14) = [-1 2 -1 2; 1 -2 1 3];
    ssa(:,:,15) = [-2 1 1 2; 2 -1 -1 3];
    ssa(:,:,16) = [-1 -1 2 2; 1 1 -2 3];
    ssa(:,:,17) = [1 -2 1 2; -1 2 -1 3];
    ssa(:,:,18) = [2 -1 -1 2; -2 1 1 3];
    % extension twin
    if strcmpi(str2,'twin')
        ntwin = 6;
        ssGroup = [ssGroup,ssGroup(end)+6];
        ssa(:,:,19) = [1 0 -1 2; -1 0 1 1];
        ssa(:,:,20) = [0 1 -1 2; 0 -1 1 1];
        ssa(:,:,21) = [-1 1 0 2; 1 -1 0 1];
        ssa(:,:,22) = [-1 0 1 2; 1 0 -1 1];
        ssa(:,:,23) = [0 -1 1 2; 0 1 -1 1];
        ssa(:,:,24) = [1 -1 0 2; -1 1 0 1];
    end
    if strcmpi(str2,'ctwin')
        % contraction twin
        ntwin = 12;
        ssGroup = [ssGroup,ssGroup(end)+6,ssGroup(end)+12];
        ssa(:,:,19) = [1 0 -1 2; -1 0 1 1];
        ssa(:,:,20) = [0 1 -1 2; 0 -1 1 1];
        ssa(:,:,21) = [-1 1 0 2; 1 -1 0 1];
        ssa(:,:,22) = [-1 0 1 2; 1 0 -1 1];
        ssa(:,:,23) = [0 -1 1 2; 0 1 -1 1];
        ssa(:,:,24) = [1 -1 0 2; -1 1 0 1];
        
        ssa(:,:,25) = [1 0 -1 1; 1 0 -1 -2];
        ssa(:,:,26) = [0 1 -1 1; 0 1 -1 -2];
        ssa(:,:,27) = [-1 1 0 1; -1 1 0 -2];
        ssa(:,:,28) = [-1 0 1 1; -1 0 1 -2];
        ssa(:,:,29) = [0 -1 1 1; 0 -1 1 -2];
        ssa(:,:,30) = [1 -1 0 1; 1 -1 0 -2];
    end

elseif any(strcmpi(str1,{'Al','Aluminum','Fcc'}))
    c_a = 1;
    nss = 12;
    ntwin = 0;
    ssGroup = [3,6,9,12];
    ssa(:,:,1) = [1 1 1; 1 0 -1];
    ssa(:,:,2) = [1 1 1; -1 1 0];
    ssa(:,:,3) = [1 1 1; 0 -1 1];
    
    ssa(:,:,4) = [-1 1 1; 0 1 -1];
    ssa(:,:,5) = [-1 1 1; -1 -1 0];
    ssa(:,:,6) = [-1 1 1; 1 0 1];
    
    ssa(:,:,7) = [-1 -1 1; -1 0 -1];
    ssa(:,:,8) = [-1 -1 1; 1 -1 0];
    ssa(:,:,9) = [-1 -1 1; 0 1 1];
    
    ssa(:,:,10) = [1 -1 1; 0 -1 -1];
    ssa(:,:,11) = [1 -1 1; 1 1 0];
    ssa(:,:,12) = [1 -1 1; -1 0 1];
    % add twin system here --------------------------------
    if strcmpi(str2,'twin')
        ntwin;
        ssGroup;
%         ssa(:,:,19) = [1 0 -1 2; -1 0 1 1];
%         ssa(:,:,20) = [0 1 -1 2; 0 -1 1 1];
%         ssa(:,:,21) = [-1 1 0 2; 1 -1 0 1];
%         ssa(:,:,22) = [-1 0 1 2; 1 0 -1 1];
%         ssa(:,:,23) = [0 -1 1 2; 0 1 -1 1];
%         ssa(:,:,24) = [1 -1 0 2; -1 1 0 1];
    end
else
    disp('no ss data for this material.');
    ssa = [];
    c_a = [];
end