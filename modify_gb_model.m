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
S = evalin('base','S;');
stepSize = evalin('base','stepSize;');

% select a point, which is closest to the position you click
h_temp = impoint(gca);
pos = h_temp.getPosition;
delete(h_temp);
[~,ind] = min(pdist2(pt_pos,pos));

delete_to_pos = [0, 0]; % set deleted point to this position

switch method
    case {'remove_point',1}
        disp('remove_point');
        
        % Set its [pt_pos] to delete_to_pos 
        pt_pos(ind,:) = delete_to_pos;
        
        % affected points
        pts_for_this_boundary = [];
        
        % from its pt_s_gb, do ...
        for ii=1:length(pt_s_gb{ind})
           igb = pt_s_gb{ind}(ii);
           
           pts_for_this_boundary = [pts_for_this_boundary, gb_s_pt{igb}];
           
           idx = find(gb_s_pt{igb}==ind);
           % Delete it from the [gb_s_pt]
           gb_s_pt{igb}(idx) = [];
           % Delete it from [H]
           H{igb}(idx) = [];
        end
 
        % Have to remove [callback], updating only does not work.
        removeNewPositionCallback(h{ind},S{ind});
        
        % Delete its [pt_s_gb], so it no longer has a gb
        pt_s_gb{ind} = [];

        % move [h]
        h{ind}.setPosition(delete_to_pos);
        % Also need to modify h,H,hline,L,G,V
        
        % Update all related points.   
        % Need to modify [L,G,V,callback,S], where [L{ind},G{ind}] required to update callback.
        % [h,H,hline] are OK.
        for ipt = pts_for_this_boundary
            L{ipt} = hline(pt_s_gb{ipt});
            G{ipt} = H(pt_s_gb{ipt});
            V{ipt} = gb_dir(pt_s_gb{ipt});
            S{ipt} = addNewPositionCallback(h{ipt}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,stepSize) , L{ipt}, G{ipt}, V{ipt}) );
        end
        
        
    case {'add_boundary',2}
        disp('add_boundary');

        
        
    case {'vertical','horizontal'}
        disp('change_gb_dir');
        % for each of its [pt_s_gb], change the [gb_dir]
        for ii=1:length(pt_s_gb{ind})
            igb = pt_s_gb{ind}(ii);
            gb_dir{igb} = method;            
        end
        
        % Update all related points.   
        % Need to modify [L,G,V,callback,S], where [L{ind},G{ind}] required to update callback.
        % [h,H,hline] are OK.
        for ii=1:length(pt_s_gb{ind})
            igb = pt_s_gb{ind}(ii)
            for jj = 1:length(gb_s_pt{igb})
                ipt = gb_s_pt{igb}(jj)
                L{ipt} = hline(pt_s_gb{ipt});
                G{ipt} = H(pt_s_gb{ipt});
                V{ipt} = gb_dir(pt_s_gb{ipt});
                S{ipt} = addNewPositionCallback(h{ipt}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,stepSize) , L{ipt}, G{ipt}, V{ipt}) );
            end
        end
        
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
assignin('base','S',S);

end