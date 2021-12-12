% load_settings(path_of_setting_file, variable_names(opt) )
% load all variables in the setting files, or selected variables.
%
% chenzhe, 2017-05-30

function [] = load_settings(pf, varargin)

d = load(pf);

if isempty(varargin)
    names = fieldnames(d);
    for ii = 1:length(names)
        assignin('base',names{ii},d.(names{ii}));
    end
else
    for ii = 1:length(varargin)
        assignin('base',varargin{ii},d.(varargin{ii}));
    end
end

end