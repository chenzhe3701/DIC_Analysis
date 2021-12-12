% chenzhe, 2018-02-09
% find ID from a map
%
% chenzhe, 2018-03-28. add [suby, subx] as output, which can be used as
% [indr, indc] for other purposes.

function [ids, indr, indcx] = find_ID_on_map(X,Y,ID,f,a)
    set(f,'currentaxes',a);
    [x,y] = getpts;
    [nr,nc] = size(ID);
    for ii = 1:size(x,1)
       [~,subx] = min(abs(X(1,:)-x(ii)));
       [~,suby] = min(abs(Y(:,1)-y(ii)));
       ids(ii) = ID(suby,subx);
       indr(ii) =  suby;
       indc(ii) =  subx;
    end
end