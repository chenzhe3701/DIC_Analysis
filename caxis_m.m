% adjust caxis for map plotted by myplotm
%
% chenzhe, 2018-04-11

function c = caxis_m(clim)

a = clim(1);
b = clim(2);

caxis([a,b]);
c = colorbar;
clim = caxis;
c_step = (clim(2)-clim(1))/size(colormap,1);
set(c,'Limits',clim+[c_step,-c_step]);

end