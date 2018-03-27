% this script is used after determine_nClusters, to look at its
% distribution
%
% chenzhe, 2018-03-25

id = zeros(1,1417);
% gD = zeros(1,1417);
nCluster = zeros(1,1417);
nClusterNew = zeros(1,1417);
for ii=1:1417
    id(ii) = struNC(ii).gID;
%     gD(ii) = gDiameter(gID==id(ii));
    nCluster(ii) = struNC(ii).nCluster;
    nClusterNew(ii) = struNC_new(ii).nCluster;
end

figure; plot(id, nClusterNew-nCluster, 'o');