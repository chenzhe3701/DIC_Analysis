
function out = color_sc(exx,cmin,cmax)

cFrom = linspace(cmin,cmax,64);
cTo = linspace(0,1,64);
out = interp1(cFrom,cTo,exx);