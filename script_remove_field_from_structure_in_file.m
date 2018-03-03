
% flexible code to remove a field from a structure in a file
%
% chenzhe, 2018-03-03

saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

for iE=2:5
    fName_c2t_result = ['WE43_T6_C1_s',num2str(iE),'_cluster_to_twin_result.mat'];
    stru = load([saveDataPath,fName_c2t_result],'stru');   
    stru = stru.stru;
    
    stru = rmfield(stru,'vInc');
    save([saveDataPath,fName_c2t_result],'stru','-append');    
end

    