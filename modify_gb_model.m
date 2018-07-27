% codes to moodify_gb_model
%
% chenzhe, 2018-05-17
%
% Things needed:
% (1) remove point
% (2) add boundary
% (3) change gb_dir

function modify_gb_model(method,varargin)

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

delete_to_pos = [nan, nan]; % set deleted point to this position

% find out the grain boundaries that are actually plotted in this AOI.
try
    all_gb_ind = evalin('base','all_gb_ind;');
catch
    all_gb_ind = 1:length(gb_s_pt);
end

switch method
    case {'check id', 0}
        disp(['point ind: ', num2str(ind)]);
        
        % from its pt_s_gb, do ...
        for ii=1:length(pt_s_gb{ind})
           disp(['boundary: ' num2str(pt_s_gb{ind}(ii))]);
        end
        
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
            
            pt_s_gb_in_aoi = intersect(pt_s_gb{ipt}, all_gb_ind);
            
            L{ipt} = hline(pt_s_gb_in_aoi);
            G{ipt} = H(pt_s_gb_in_aoi);
            V{ipt} = gb_dir(pt_s_gb_in_aoi);
            
            removeNewPositionCallback(h{ipt},S{ipt});   % need to remove first, otherwise it keep updating
            
            S{ipt} = addNewPositionCallback(h{ipt}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,p,ipt,stepSize) , L{ipt}, G{ipt}, V{ipt}) );
        end
        
        
    case {'add_boundary',2}
        try
            npts_add = varargin{1};
        catch
            npts_add = 3;
        end
        disp('add_boundary');
        % make addition to [gb_dir, gb_s_pt, pt_pos, pt_s_gb]
        ngb = length(gb_dir);
        npt = size(pt_pos,1);
        gb_dir{ngb+1} = 'horizontal';
        gb_s_pt{ngb+1} = npt + [1:npts_add];
        pt_pos(npt+[1:npts_add],:) = pos + [stepSize,stepSize].*[linspace(0,50*npts_add,npts_add)]';
        for ii=1:npts_add
            pt_s_gb{npt+ii} = ngb+1;
        end
        
        % make addition to these 
        % (1) plot all the control points --> hangle: h{i}
        for ii = npt+1 : npt+npts_add
            h{ii} = impoint(gca, pt_pos(ii,1), pt_pos(ii,2));
            setColor(h{ii},'r');
        end
        % (2) for each boundary, group all its impoint handles --> gb_s_pt_group{j} = H{j} = {h{j1}, h{j2}, ... }
        for jj = ngb+1
            H{jj} = h(gb_s_pt{jj}); % or, looks like this is the same { h{ gb_s_pt{jj} } }
        end
        % (3) for each boundary, plot the line and record the handle: hline{j}
        for jj = ngb+1
            hline{jj} = plot_spline_line(H{jj}, gb_dir{jj}, stepSize);
        end
        
        % (4) for each point, group its related grain boundaries pt_s_gb_group{i} = L{i} = {hline{i1}, hline{i2}, ...}
        % group this point's grain boundaries' points, pt_s_gb_s_pt_group{i} = G{i} = {H{i1}, H{i2}, ...}
        % gropu this point's grain boundaries' direction, pt_s_gb_sdir{i} = V{i} = {gb_dir{i1}, gb_dir{i2}, ...}
        % addNewPositionCallback with L,G,V
        for ii = npt+1 : npt+npts_add
            L{ii} = hline(pt_s_gb{ii});
            G{ii} = H(pt_s_gb{ii});
            V{ii} = gb_dir(pt_s_gb{ii});
            S{ii} = addNewPositionCallback(h{ii}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,p,ii,stepSize) , L{ii}, G{ii}, V{ii}) );
        end

        try
           all_pts_ind = evalin('base','all_pts_ind;');
           all_pts_ind = [all_pts_ind, npt+1 : npt+npts_add]; 
           assignin('base','all_pts_ind',all_pts_ind);
        end
        try
            all_gb_ind = evalin('base','all_gb_ind;');
            all_gb_ind = [ all_gb_ind, ngb+1];
            assignin('base','all_gb_ind',all_gb_ind);
        end
        
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
            igb = pt_s_gb{ind}(ii);
            for jj = 1:length(gb_s_pt{igb})
                ipt = gb_s_pt{igb}(jj);
                
                pt_s_gb_in_aoi = intersect(pt_s_gb{ipt}, all_gb_ind);
                
                L{ipt} = hline(pt_s_gb_in_aoi);
                G{ipt} = H(pt_s_gb_in_aoi);
                V{ipt} = gb_dir(pt_s_gb_in_aoi);
                
                removeNewPositionCallback(h{ipt},S{ipt});   % need to remove first, otherwise it keep updating
                
                S{ipt} = addNewPositionCallback(h{ipt}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,p,ipt,stepSize) , L{ipt}, G{ipt}, V{ipt}) );
            end
        end
    case{'merge',4}
        
        % select a second point, which is closest to the position you click
        h_temp = impoint(gca);
        pos = h_temp.getPosition;
        delete(h_temp);
        [~,ind_merge_to] = min(pdist2(pt_pos,pos));
        
        disp('remove_point');
        
        % Set its [pt_pos] to delete_to_pos 
        pt_pos(ind,:) = delete_to_pos;
        
        % affected points
        pts_for_this_boundary = [];
        
        % from its pt_s_gb, do ...
        for ii=1:length(pt_s_gb{ind})
           igb = pt_s_gb{ind}(ii);
           
           % add grain boundary to pt_s_gb{ind_merge_to}
           pt_s_gb{ind_merge_to} = [pt_s_gb{ind_merge_to}, igb];
           
           idx = find(gb_s_pt{igb}==ind);
           % From the [gb_s_pt], change ind to ind_merge_to
           gb_s_pt{igb}(idx) = ind_merge_to;
           % Modify the gb's handle group of impoints
           H{igb}(idx) = h(ind_merge_to);
           
           pts_for_this_boundary = [pts_for_this_boundary, gb_s_pt{igb}];
           
        end

        % Have to remove [callback], updating only does not work.
        removeNewPositionCallback(h{ind},S{ind});
        
        % Delete its [pt_s_gb], so it no longer has a gb
        pt_s_gb{ind} = [];


        % Also need to modify h,H,hline,L,G,V
        
        % Update all related points.   
        % Need to modify [L,G,V,callback,S], where [L{ind},G{ind}] required to update callback.
        % [h,H,hline] are OK.
        for ipt = pts_for_this_boundary
            
            pt_s_gb_in_aoi = intersect(pt_s_gb{ipt}, all_gb_ind);
            
            L{ipt} = hline(pt_s_gb_in_aoi);
            G{ipt} = H(pt_s_gb_in_aoi);
            V{ipt} = gb_dir(pt_s_gb_in_aoi);
            
            removeNewPositionCallback(h{ipt},S{ipt});   % need to remove first, otherwise it keep updating
            
            S{ipt} = addNewPositionCallback(h{ipt}, @(p) cellfun(@(x,y,z) update_spline_line_hv(x,y,z,p,ipt,stepSize) , L{ipt}, G{ipt}, V{ipt}) );
        end
        
        % move [h] to deleted position
        h{ind}.setPosition(delete_to_pos);
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