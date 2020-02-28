function ms_out = MS_remove_data_sandbox(cfg_in,ms_in, cells_to_remove)
%% MS_remove_data
%  MS_remove_data will take in a miniscope data structure 'ms_in' and
%  remove the specified cells while performing some checks to ensure that everything lines up.
%
%   Example:
%       ms_out = MS_remove_data(ms_in, cells_to_remove)
%           The ms data coming in should be in the form of segmented
%          miniscope data (suggest using 'MS_segment_ms') and will remove
%          all cells listed in cells_to_remove which is an array of indices.
%
%   Inputs:
%       - ms_in: [struct] the ms data structure from
%       msGenerateVideoObj_ConcatenateBatch*. Should be pre-segmented into
%       cells per segment.
%       - cells_to_remove: [1xN] an array of indicies to remove. ex:
%       [1,4,19]
%
%   Outputs:
%       - ms_out: [struct] with removed data segment cells.
%
%
%  * in https://github.com/torilynntemple/MiniscopeAnalysisConcatenation by
%   Tori-Lynn Temple (Williams lab at McGill)
%
%   EC 2019-12-16
%       - initial sandbox with basic functions
%% inialize

cfg_def = [];
cfg_def.user_fields = {}; % user fields. 

cfg = ProcessConfig(cfg_def, cfg_in); 

if ~iscell(ms_in.RawTraces)
    error('ms_in data is not in cells.  Probably has not been segmented.  See MS_segment_ms_data')
end


ms_out = ms_in; % copy to maintain all other contents

% just for later check
segments_in = 1:length(ms_in.RawTraces);
keep_idx = segments_in(~ismember(segments_in, cells_to_remove));


%% remove all cells in the data fields with pre-segmented data
% fields_to_alter = {'time', 'RawTraces', 'FiltTraces', 'frameNum', 'vidNum'};

known_cell_num = size(ms_in.RawTraces,1); % should always be the correct number of cells for the number of segments.
fields = fieldnames(ms_in);
for iF = 1:length(fields)
    field_size = size(ms_in.(fields{iF}));
    
    cell_idx = find(field_size == known_cell_num,1);
    
    
    if ~isempty(cell_idx)
        
        for iC = sort(cells_to_remove, 'descend') % go backards or else everything is off and will not remove the right cells.
            ms_out.(fields{iF})(iC) = [];
            
            
            fprintf('Removing cell: %0.f in %s\n', iC, fields{iF});
        end
    end
end



%% check in the inputs vs the outputs

% compare how many cells remain vs how many should be removed. 
if length(ms_out.RawTraces) ~= length(keep_idx)
    error('Size of ms_in does not match the size of ms_out given cells_to remove')
end

% check the lengths align
for iK = length(ms_out.RawTraces):-1:1
    checks(iK) = isequal(ms_in.RawTraces{keep_idx(iK)} , ms_out.RawTraces{iK});
end

%if the lengths are off, panic. 
if sum(checks) ~= length(checks)
    error('Size of ms_in does not match the size of ms_out given cells_to remove')
end


%% clean up
ms_out.format = [ms_out.format sprintf('cell %.0f removed', cells_to_remove)];

if isfield(ms_out, 'history')
    ms_out.history.fun_name{end+1} = {sprintf('cell %.0f removed', cells_to_remove)};
    ms_out.history.cfg{end+1} = cfg;
else
    ms_out.history.fun_name = {sprintf('cell %.0f removed', cells_to_remove)};
    ms_out.history.cfg = {cfg};
end

