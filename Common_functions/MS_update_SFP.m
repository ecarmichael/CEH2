function ms_out = MS_update_SFP(cfg_in, ms_in)
%% MS_update_SFP: goes through the binary v
%
%
%
%    Inputs: 
%     -
%
%
%
%    Outputs: 
%     -
%
%
%
%
% EC 2020-03-04   initial version 
%
%
%% initialize

cfg_def = [];
cfg_def.field_to_use = 'Binary'; % which ms_in subfield to use.
cfg_def.fnc = '=='; % which function to use. '==', '>', '>=', '<', '<='
cfg_def.remove_val = 0; % value to use in selection. 0 by default. 
cfg_def.replace_val = 0; % what do you want to replace it with? 
cfg_def.centroid_val = NaN; 
cfg = ProcessConfig2(cfg_def, cfg_in); 


% if ~isfield('ms_in', 'binary')
%     error('requires binary field in ms_in')
% end

ms_out = ms_in; 
%% cycle through cells and determine if it is active during this period.  If not replace with 0
removed_cells = []; 

for iC = 1:size(ms_in.(cfg.field_to_use),2)
   
    if sum(ms_in.(cfg.field_to_use)(:,iC)) == 0
        ms_out.SFPs(:,:,iC) = cfg.replace_val; 
%         ms_out.Centroids(iC,:) = cfg.centroid_val; 
        
        removed_cells = [removed_cells, iC]; % used for tracking later. 
    end
    

end



%% clean up

fprintf('<strong>%s</strong>: %d cells reaplced with %d from ms.SFP due to inactivity\n', mfilename, length(removed_cells), cfg.replace_val); 


