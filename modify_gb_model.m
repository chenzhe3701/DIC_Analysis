% codes to moodify_gb_model
%
% chenzhe, 2018-05-17
%
% Things needed:
% (1) remove point
% (2) add boundary
% (3) change gb_dir

function modify_gb_model(method)

gb_dir = evalin('base','gb_dir;');
gb_s_pt = evalin('base','gb_s_pt;');
pt_pos = evalin('base','pt_pos;');
pt_s_gb = evalin('base','pt_s_gb;');
tripleLookup = evalin('base','tripleLookup;');
h = evalin('base','h;');
H = evalin('base','H;');
hline = evalin('base','hline;');
L = evalin('base','L;');
G = evalin('base','G;');
V = evalin('base','V;');

% select a point, which is closest to the position you click
h_temp = impoint(gca);
pos = h_temp.getPosition;
delete(h_temp);
[~,ind] = min(pdist2(pt_pos,pos));

switch method
    case {'remove_point',1}
        disp('remove_pont');
        % remove point.      
        % By setting its [pt_pos] to [-100,-100] 
        pt_pos(ind,:) = [-100 -100];
        % from its pt_s_gb, delete it from the [gb_s_pt]
        for ii=1:length(pt_s_gb{ind})
           igb = pt_s_gb{ind}(ii);
           gb_s_pt{igb}(gb_s_pt{igb}==ind) = [];
        end
        % delete its [pt_s_gb], so it no longer has a gb
        pt_s_gb{ind} = [];
        
    case {'add_boundary',2}
        disp('add_boundary');
        
    case {'change_gb_dir',3}
        disp('change_gb_dir');
end



assignin('base','gb_dir',gb_dir);
assignin('base','gb_s_pt',gb_s_pt);
assignin('base','pt_pos',pt_pos);
assignin('base','pt_s_gb',pt_s_gb);
assignin('base','tripleLookup',tripleLookup);
assignin('base','h',h);
assignin('base','H',H);
assignin('base','hline',hline);
assignin('base','L',L);
assignin('base','G',G);
assignin('base','V',V);

end