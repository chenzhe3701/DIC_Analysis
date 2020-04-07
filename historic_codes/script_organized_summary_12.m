
% based on 11, temp.
% Comment, 2019-08-16:
% Looks like this only consider if neighbor grain has a twin, but do not consider which grain boundary the twin intersects.   
% When look at mPrime, it simply choose the basal in neighbor with the highest SF.   
% So, looks like this is not a very good code.  The code "show_manual_label_..." should be a better code. 

iE = 2;
variableNames3 = {'iE','ID','iTwin','ID_neighbor','basal_SF_nb','mPrime','twinsNbTF'};
T3 = cell2table(cell(0,length(variableNames3)));
T3.Properties.VariableNames = variableNames3;

for iS = 1:length(struCell{iE})
    ID_current = struCell{iE}(iS).gID
    
    % [*] Calculate the max basal Schmid factor of this grain
    ind = find(ID_current==gID);
    euler_current = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
%     [abs_schmid_factor, ~, burgersXY] = trace_analysis_TiMgAl(euler_current, [0 0 0], [0 0 0], stressTensor, 'Mg', 'twin');
%     basal_SF = max(abs_schmid_factor(1:3,2));
%     twin_SF = max(abs_schmid_factor(19:24,2));
    
    
    % All neighboring grains
    nNeighbors = gNNeighbors(ind);
    ID_neighbors = gNeighbors(ind, 1:nNeighbors);
    
    
    % find out new activeTS
    activeTS = logical(sum(struCell{iE}(iS).cTrueTwin,1));
    

    % % %
    % find all twinned grains
    all_twinned_grains = unique(tb_gNum);
    all_involved_grains = unique([floor(tb_gbNum/10000); mod(tb_gbNum,10000)]);
    
    
    % If this grain is twinned.  Then, for each active
    if any(activeTS)
        for iTwin = 1:length(activeTS)
            tgbs = struCell{iE}(iS).tGb{iTwin};    % find all of this twin variant's grain boundaries   
            
            % Find tgbs that first appear
            for ie = iE-1:-1:2
                tgbs_earlier = struCell{ie}(iS).tGb{iTwin};
                if ~isempty(tgbs_earlier)
                    tgbs = tgbs_earlier;
                end
            end

            % Method-1
            if ~isempty(tgbs)
                % This mean this twin variant has intersecting grain boundary.  
                
                % Find out this twin's neighbors
                nTwinsNbs = length(tgbs);   
                twinsNbs = [];
                for iTwinsNb = 1:nTwinsNbs
                    tgb_current = tgbs(iTwinsNb);
                    pairs = [floor(tgb_current/10000), mod(tgb_current,10000)];
                    ID_neighbor = pairs(~ismember(pairs,ID_current));   % find neighbor ID
                    if length(ID_neighbor)==1
                    twinsNbs(iTwinsNb) = ID_neighbor;                   % record this twins neighbor.  Note iE=2, ID_442, twinGB does not have current ID, which is wrong and have to debug.
                    end
                end
                nonTwinsNbs = ID_neighbors(~ismember(ID_neighbors,twinsNbs));   % neighbors that are not related to this twin
                
                % Loop for all grain neighbors
                for iNb = 1:length(ID_neighbors)
                    ID_neighbor = ID_neighbors(iNb);
                    ind = find(ID_neighbor==gID);
                    euler_neighbor = [gPhi1(ind),gPhi(ind),gPhi2(ind)];
                    [schmidFactorG1, schmidFactorG2, mPrimeMatrix, resBurgersMatrix, mPrimeMatrixAbs, resBurgersMatrixAbs] = calculate_mPrime_and_resB(euler_current, euler_neighbor, stressTensor, [1 0 0], sampleMaterial, 'twin');
                    
                    % if the grain neighbor is a twin neighbor
                    if ismember(ID_neighbor,twinsNbs)
                        % Check if this twin neighbor is basal slipped / twinned
                        twinsNbTF = true;
                        if ~ismember(ID_neighbor,all_twinned_grains)
                            [basal_SF_nb, iBasal] = max(schmidFactorG2(1:3));      % max basal SF in the neighboring grain
                            mPrime = mPrimeMatrixAbs(19:24,1:3);
                            mPrime = mPrime(iTwin,iBasal);      % max of this twin w.r.t a basal slip in neighbor
                            
                            T3 = [T3; {iE, ID_current, iTwin, ID_neighbor,basal_SF_nb, mPrime, twinsNbTF}];
                        else

                        end
                    else
                        % Check if this non-twin-neighbor is basal slipped / twinned  
                        twinsNbTF = false;
                        if ~ismember(ID_neighbor,all_twinned_grains)
                            [basal_SF_nb, iBasal] = max(schmidFactorG2(1:3));      % max basal SF in the neighboring grain
                            mPrime = mPrimeMatrixAbs(19:24,1:3);
                            mPrime = mPrime(iTwin,iBasal);      % max of this twin w.r.t a basal slip in neighbor
                            
                            T3 = [T3; {iE, ID_current, iTwin, ID_neighbor,basal_SF_nb, mPrime, twinsNbTF}];
                        else
                            
                        end
                    end
                end
                
                
            else
                disp(['twinned, but do not have intersecting grain boundary:',newline,'iE:',num2str(iE),', iS:',num2str(iS),', ID:',num2str(ID_current),', TS:',num2str(iTwin)]);
            end
            
        end
    else
        % Not twinned
    end
    
    
end


%%  Summarize basal_SF_Nb vs. mPrime, between 2 groups: (1) neighbor twinned, (2) neighbor not twinned.
%  The difference is not big.  And high basal tend to have high m', but what does it mean?
close all
inds = T3.twinsNbTF==1;
T4 = T3(inds,:)
figure;hist(T3.basal_SF_nb)
figure;hist(T4.basal_SF_nb)
figure;hist(T3.mPrime)
figure;hist(T4.mPrime)
figure;plot(T3.basal_SF_nb,T3.mPrime,'.')
figure;plot(T4.basal_SF_nb,T4.mPrime,'.')



