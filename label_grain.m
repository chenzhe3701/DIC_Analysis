% chenzhe, 2018-03-30
% 
% make a label on map with the ID interested in, to help locate the grain.
%
% hf = label_map_with_ID(X,Y,ID,hf)
% X, Y, ID are matrice of the same size: coordinates and ID_map.
% hf is the handle_of_figure

function hf = label_grain(target_id,X,Y,ID,hf)

set(0,'currentfigure',hf);

hold on;

ID_reduced = ID(1:10:end,1:10:end);
X_reduced = X(1:10:end,1:10:end);
Y_reduced = Y(1:10:end,1:10:end);
z = gca;
z = z.ZLim(2);
for ii = 1:length(target_id)
    id = target_id(ii);
    ind = (ID_reduced==id);
    x = mean(X_reduced(ind));
    y = mean(Y_reduced(ind));
    text(x,y,z,num2str(id));
end

hold off;

end