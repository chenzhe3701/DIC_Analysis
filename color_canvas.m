% [RGB,R,G,B] = color_canvas(nR,nC, mapNumber)
% output a matrix of size [nR,nC], that has a color gradient
% can define different kinds of maps using mapNumber
% chenzhe, 2018-07-24

function [RGB,R,G,B] = color_canvas(nR,nC, mapNumber)

if ~exist('mapNumber','var')
    mapNumber = 1;
end

switch mapNumber
    case 1
        % [R K; G B]
        R1 = [1 0; 0 0];
        G1 = [0 0; 1 0];
        B1 = [0 0; 0 1];
    case 2
        
        % [R M; G B]
        R1 = [1 1; 0 0];
        G1 = [0 0; 1 0];
        B1 = [0 1; 0 1];
end

% CMYK
% R1 = [0 1; 1 0];
% G1 = [1 0; 1 0];
% B1 = [1 1; 0 0];


[x1,y1] = meshgrid(0:1,0:1);
[x2,y2] = meshgrid(linspace(0,1,nC),linspace(0,1,nR));

R = interp_data(x1,y1,R1,x2,y2,[],'interp','linear');
G = interp_data(x1,y1,G1,x2,y2,[],'interp','linear');
B = interp_data(x1,y1,B1,x2,y2,[],'interp','linear');

RGB = cat(3,R,G,B);

% if ~exist('n','var')
%     n = 65535;
% end
% 
% [IND,map] = rgb2ind(RGB,n);