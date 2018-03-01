% chenzhe, 2018-02-09
% find ID from a map

function tNote = disable_twin(X,Y,ID,clusterNumMap,tNote,f,a)
    set(f,'currentaxes',a);
    [x,y] = getpts;
    [nr,nc] = size(ID);
    for ii = 1:size(x,1)
       [~,subx] = min(abs(X(1,:)-x(ii)));
       [~,suby] = min(abs(Y(:,1)-y(ii)));
       idNum(ii) = ID(suby,subx);
       cNum(ii) = clusterNumMap(suby,subx);
    end
    tNote.disable = [tNote.disable;[idNum(:),cNum(:)]];
    tNote.disable = unique(tNote.disable,'rows');
end