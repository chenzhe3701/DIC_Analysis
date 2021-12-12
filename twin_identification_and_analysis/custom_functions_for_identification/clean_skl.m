% chenzhe, 2019-02-04, custom function to remove small branches < minLength
% from skl (logical matrix)
%
% algorithm:
% (1) find out skl=skeleton, B = branch points, E = end points
% (2) for each end point, find its distance to closest branch point
% if dist < th, then remove all skl points with dist<th to this End point
% B = original branch points for record, maybe useful

function [skl, B_old] = clean_skl(skl,minLength)

skl = logical(skl);

B =  bwmorph(skl,'branchpoints');   % branch points
B_old = B;                          % keep original branch points
E  = bwmorph(skl,'endpoints');      % end points
inds = find(E);
D = bwdistgeodesic(skl,B);  % D = bwdistgeodesic(skl,false(size(B)));  % for debug
d_from_ends = D(inds);
[val,ind] = min(d_from_ends);   % find the smallest distance & ind of the End point to a Branch point

count = 0;
while (val<minLength)%&&(count < 40)
    d =  bwdistgeodesic(skl, inds(ind));    % find the distance of all SKL points to this End point
    toRemove = d<val;
    skl(toRemove) = 0;  % myplot(skl);
    
    B =  bwmorph(skl,'branchpoints');   % branch points
    E  = bwmorph(skl,'endpoints');      % end points
    inds = find(E);
    D = bwdistgeodesic(skl,B);  % D = bwdistgeodesic(skl,false(size(B)));  % for debug
    d_from_ends = D(inds);
    [val,ind] = min(d_from_ends);
    
    count = count + 1;
end

%% Here is the code to clean, use while to iteratively remove, so some branches from end is no longer a branch after removing one branch.  But this will also cause troubles...
if 0
    close all;
    minLength = round(min(size(trueTwinMapL)) * 0.05)
    skl = bwskel(imbinarize(trueTwinMapL));
    myplot(skl);
    
    B =  bwmorph(skl,'branchpoints');   % branch points
    B_old = B;                          % keep original branch points
    E  = bwmorph(skl,'endpoints');      % end points
    inds = find(E);
    D = bwdistgeodesic(skl,B);  % D = bwdistgeodesic(skl,false(size(B)));  % for debug
    d_from_ends = D(inds);
    [val,ind] = min(d_from_ends);   % find the smallest distance & ind of the End point to a Branch point
    
    count = 0;
    while (val<minLength)%&&(count < 40)
        d =  bwdistgeodesic(skl, inds(ind));    % find the distance of all SKL points to this End point
        toRemove = d<val;
        skl(toRemove) = 0;  % myplot(skl);
        
        B =  bwmorph(skl,'branchpoints');   % branch points
        E  = bwmorph(skl,'endpoints');      % end points
        inds = find(E);
        D = bwdistgeodesic(skl,B);  % D = bwdistgeodesic(skl,false(size(B)));  % for debug
        d_from_ends = D(inds);
        [val,ind] = min(d_from_ends);
        
        count = count + 1;
    end
    
    myplot(skl);
end

%% faster, intial version
if 0
    close all;
    minLength = round(min(size(trueTwinMapL)) * 0.05)
    skl = bwskel(imbinarize(trueTwinMapL));
    myplot(skl);
    
    B =  bwmorph(skl,'branchpoints');   % branch points
    B_old = B;                          % keep original branch points
    E  = bwmorph(skl,'endpoints');      % end points
    inds = find(E);
    
    D = bwdistgeodesic(skl,B);  % D = bwdistgeodesic(skl,false(size(B)));  % for debug
    d_from_ends = D(inds);
    toRemove = false(size(skl));
    
    for ii = 1:length(inds)
        this_dist = d_from_ends(ii);
        if this_dist < minLength
           d = bwdistgeodesic(skl, inds(ii));
           toRemove = (d<this_dist)|toRemove;
        end
    end
    skl(toRemove) = 0;
    
    myplot(skl);
end



end