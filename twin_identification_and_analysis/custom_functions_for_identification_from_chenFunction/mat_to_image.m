% out = mat_to_image(M_in, m_threshold, setting)
% chenzhe, 2018-03-31%
%
% based on color_sc, which is a code I wrote early March but did not have
% time to document.
%
% This is similar to built-in mat2gray(). However, you can specify the
% range of the input data that you'd like to interp to the output limits.
%
% Convert input data matrix M within [m_low_th, m_high_th] into data ranged
% in, such as [0 1], [0 255] etc, linearly, based on the setting selected,
% so that it can be used as index in format such as RGB.
%
% Input data at [m_low_th, m_high_th] correspond to output data at [0 1],
% or [0 255] etc.

function out = mat_to_image(M_in, m_threshold, setting)

m_low_th = m_threshold(1);
m_high_th = m_threshold(2);

switch setting
    case {'index'}
        out_low = 0;
        out_high = 1;
        resolution = 2;
        
        cFrom = linspace(m_low_th, m_high_th, resolution);
        cTo = linspace(out_low, out_high, resolution);
        out = double( interp1(cFrom, cTo, M_in));
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

