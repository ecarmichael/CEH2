function ms_data_out = MS_characterize_trace(cfg_in,ms_data_in)
%% MS_Ca_trace_stats: computes a variety of characterists for each calcium trace
%   - stats.(cfg.data_type).var = variance
%   - stats.(cfg.data_type).z_max = maxium of the zscore signal. 
%   - stats.(cfg.data_type).min_max = raw min and max [min max]; 
%
%
%    Inputs: 
%     - cfg_in: [struct] user configurations 
%           - cfg.data_type = 'RawTraces' could also be 'BinaryTraces' or 'Filtertraces'
%
%
%     - ms_data_in: [struct] the miniscope data structure following
%     preprocessing
%
%
%
%    Outputs: 
%     - ms_data_out: [struct] same format with additional fields. 
%
%    ToDo:
%       - ??   
%
%
% EC 2020-01-28   initial version 
%
%% initialize

cfg_def = [];
cfg_def.data_type = 'RawTraces'; 

cfg = ProcessConfig(cfg_def, cfg_in);

%% get some statistics. 


ms_data_out = ms_data_in; 
for iCell = size(ms_data_in.(cfg.data_type),2):-1:1
    % get the max of the zscore.  useful for finding cells with no major spikes. 
ms_data_out.stats.(cfg.data_type).z_max(iCell) = max(zscore(ms_data_in.(cfg.data_type)(:,iCell)));

% get the raw min and max
ms_data_out.stats.(cfg.data_type).min_max(:,iCell) = [min(ms_data_in.(cfg.data_type)(:,iCell)), max(ms_data_in.(cfg.data_type)(:,iCell))];

% get the variance. 
ms_data_out.stats.(cfg.data_type).var(iCell) = var(ms_data_in.(cfg.data_type)(:,iCell));

end





