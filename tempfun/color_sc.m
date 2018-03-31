% chenzhe, 2018-03-31
% based on color_sc, which is a code I wrote early March but did not have
% time to document.
%
% Convert input data matrix M into data ranged in, such as [0 1], [0 255]
% etc, linearly, so that it can be used as index in format such as RGB.
%
% Input data at [m_low_th, m_high_th] correspond to output data at [0 1],
% or [0 255] etc.
%
% Replaced with mat_to_image().

function out = color_sc(M_in, m_low_th, m_high_th, setting)

switch setting
    case {'index'}
        out_low = 0;
        out_high = 1;
        resolution = 2;
        
        cFrom = linspace(m_low_th, m_high_th, resolution);
        cTo = linspace(out_low, out_high, resolution);
        out = interp1(cFrom, cTo, M_in);
    case {'uint8'}
        out_low = 0;
        out_high = 255;
        resolution = 2;
        
        cFrom = linspace(m_low_th, m_high_th, resolution);
        cTo = linspace(out_low, out_high, resolution);
        out = uint8( interp1(cFrom, cTo, M_in));
        
    case {'uint16'}
        out_low = 0;
        out_high = 65535;
        resolution = 2;
        
        cFrom = linspace(m_low_th, m_high_th, resolution);
        cTo = linspace(out_low, out_high, resolution);
        out = uint16( interp1(cFrom, cTo, M_in));
        
end

