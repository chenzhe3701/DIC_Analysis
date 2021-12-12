% chenzhe, 2016-3-29
% chenzhe, 2016-4-4 revised.
%
% [dist, uniqueGB] = city_block(uniqueGB, loop_limit), uniqueGB is a matrix,
% boundary is number, interior is 0 or nan
% calculate city block distance map from uniqueGB, for maximum of
% loop_limit loops.
% if no input for loop_limit, then set it = 10000 loops
% if loop_limit = -1, then don't calculate uniqueGB to make it a little bit faster
%
% output: (1) dist is the matrix containing city-block distance
% (2) uniqueGB is a matrix showing the id of the uniqueGB to which each
% pixel's city-block-distance is calculated from.
%
% chenzhe, 2018-09-18, add check, in case input contains only zero.
%
% chenzhe, 2019-01-20, I think it can be replace with built-in function.  

function [dist_man, uniqueGB] = city_block(uniqueGB,varargin)

if isempty(varargin)
    method = 'cityblock';
else
    method = varargin{1};
end
if any(uniqueGB(:))
    [dist_man, idx] = bwdist(uniqueGB,method);
    uniqueGB = uniqueGB(idx);
else
    dist_man = zeros(size(uniqueGB));
    warning('input matrix is all 0');
end

% if nansum(uniqueGB(:))==0
%     dist_man = zeros(size(uniqueGB));
%     return;
% end
% 
% uniqueGB(uniqueGB==0) = nan;
% dist_man = uniqueGB;
% dist_man(isnan(dist_man))=0;
% dist_man = double(logical(dist_man));
% [nR,nC] = size(uniqueGB);
% 
% f = [0 1 0; 1 0 1; 0 1 0];
% loop = 1;
% if isempty(varargin)
%     loop_limit = 10000;
%     fast = 0;
% elseif varargin{1} == -1
%     loop_limit = 10000;
%     fast = 1;
% else
%     loop_limit = varargin{1};
%     fast = 0;
% end
% 
% while loop<=loop_limit
%     pt_left = sum((dist_man(:)==0));
%     
%     if mod(loop, 100)==0
%         disp(['complete pct: ',num2str(1-pt_left/nR/nC)]);
%     end
%     
%     if pt_left==0
%         disp(['complete']);
%         break
%     end
%     
%     loop = loop + 1; % add 1 here, because initial dist is (0,1) instead of (nan,0)
%     if mod(loop,100)==0
%         disp(['current g.b. thickness: ',num2str(loop)]);    
%     end
%     
%     newGB = filter2(f,dist_man,'same') & (~dist_man);
%     if fast == 1
%         dist_man(newGB) = loop;     % if fast == 1, then run a little bit faster, only determine city-block distance, do not grow the uniqueGB map. 
%       else
%         [indR,indC] = find(newGB);
%         indR1 = max(indR-1,1);
%         indR2 = min(indR+1,nR);
%         indC1 = max(indC-1,1);
%         indC2 = min(indC+1,nC);
%         ind = sub2ind([nR,nC],indR,indC);
%         ind4 = sub2ind([nR,nC],[indR1;indR2;indR;indR],[indC;indC;indC1;indC2]);
%         ind4 = reshape(ind4,[],4);
%         
%         dist_man(ind) = loop;
%         uniqueGB(ind) = mode(uniqueGB(ind4),2);
%     end
%     
% end
% disp(['final g.b. thickness: ',num2str(loop)]);
% 
% dist_man = dist_man - 1;    % decrease by 1, because initial dist is (0,1) instead of (nan,0)
% uniqueGB(isnan(uniqueGB))=0;

end