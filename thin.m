% function [I, thinable] = thin(I, n)
% 
% Thin by the composite structuring element
% B1 = [1 1 1; 0 1 0; -1 -1 -1]
% B2 = [1 1 0; 1 1 -1; 0 -1 -1]
% and all their rotations, for 'n' times.
%
% If input is 3x3, thinable = if the central pixel can be thinned. 
%
% Note I have two expressions of the code which should be the same.
% But the result for thinning 'ones(3)', or maybe more generally, 
% for thinning matrices with '1's on image edge, will/might be different.
%
% Maybe Matlab have special handling of the situation.
%
% chenzhe, 2018-05-19


function [I, thinable] = thin(I, nTimes)

if ~exist('nTimes','var')
    nTimes = 1;
end

% make into 0 or 1
I(isnan(I)) = 0;
I(I~=0) = 1;

option = 1;
if option==1
    I = bwmorph(I,'thin',nTimes);
    
elseif option==2
    B1 = [1 1 1; 0 1 0; -1 -1 -1];
    B2 = [1 1 0; 1 1 -1; 0 -1 -1];
    
    I_old = zeros(size(I));
    iLoop = 0;
    
    while (iLoop < nTimes)&&(sum(I(:)-I_old(:))>0)
        for ii=1:4
            I = I& ~bwhitmiss(I, rot90(B1,ii));
            I = I& ~bwhitmiss(I, rot90(B2,ii));
        end
        I_old = I;
    end
end

% If central pixel can be thinned
[nR,nC] = size(I);
if (nR==3)&&(nC==3)
    thinable = I(5)==0;
else
    thinable = [];
end

end