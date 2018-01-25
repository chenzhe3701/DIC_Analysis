% [r,c] = find_fov(targetX,targetY)
% chenzhe, 2018-01-23
% function to find the row & column index of a fov
function [r_one_base,c_one_base] = find_fov(targetX,targetY)

[f, p] = uigetfile('D:\WE43_T6_C1_insitu_compression\stitched_img\translations_searched_vertical_stop_0.mat','find the translation file');
load([p,'\',f]);

[nR,nC] = size(transX);
rowOK = zeros(nR,nC);
colOK = zeros(nR,nC);

for iR = 1:nR
    rowOK(iR, find(transX(iR,:)<targetX,1,'last')) = 1;
end
for iC = 1:nC
   colOK(find(transY(:,iC)<targetY,1,'last'), iC) = 1; 
end
bothOK = rowOK&colOK;
[r_one_base,c_one_base] = find(bothOK)