% chenzhe, 2018-02-09
% find ID from a map

function ids = find_ID_on_map(X,Y,ID,f,a)
    set(f,'currentaxes',a);
    [x,y] = getpts;
    [nr,nc] = size(ID);
    for ii = 1:size(x,1)
       [~,subx] = min(abs(X(1,:)-x(ii)));
       [~,suby] = min(abs(Y(:,1)-y(ii)));
       ids(ii) = ID(suby,subx);
    end
end