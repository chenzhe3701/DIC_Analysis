% misorientationStruct = construct_misorientation_structure(neighborStruct, gPhi1, gPhi, gPhi2)
%
% Input neighborStruct.g1 has nn grains. 
% [gPhi1,gPhi,gPhi2] are the Euler angles for these nn grains.  The sequence should be the same.
% neighborStruct.g2 is cell array containing g1's neighbors.
%
% The output misorientationStruct is neighborStruct with addition field
% neighborStrut.misorientation{grain-1}(grain-2s)
%
% Zhe Chen, 2015-08-04 revised.
% chenzhe, 2018-07-30, modify so that if grain has no neighbors, neighborStruct.misorientation{iGrain1}(iGrain2) = nan;

function misorientationStruct = construct_misorientation_structure(neighborStruct, gPhi1, gPhi, gPhi2)

for iGrain1=1:length(neighborStruct.g1)
    phi1G1 = gPhi1(iGrain1);        % grain 1 phi's
    phiG1 = gPhi(iGrain1);
    phi2G1 = gPhi2(iGrain1);
    
    for iGrain2=1:length(neighborStruct.g2{iGrain1})
        grain2Index = find(neighborStruct.g1==neighborStruct.g2{iGrain1}(iGrain2));     % It's important to find INDEX, because grain ID numbers are not consecutive
        phi1G2 = gPhi1(grain2Index);        % grain 2 phi's
        phiG2 = gPhi(grain2Index);
        phi2G2 = gPhi2(grain2Index);
        
        if ~isempty(grain2Index)
        neighborStruct.misorientation{iGrain1}(iGrain2) = calculate_misorientation_hcp([phi1G1, phiG1, phi2G1], [phi1G2, phiG2, phi2G2]);
        else
        neighborStruct.misorientation{iGrain1}(iGrain2) = nan;
        end
    end
end
misorientationStruct = neighborStruct;