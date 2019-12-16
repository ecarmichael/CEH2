function ms_out = MS_remove_data_sandbox(ms_in, cells_to_remove)
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

if ~iscell(ms_in.RawTraces)
    error('ms_in data is not in cells.  Probably has not been segmented.  See MS_segment_ms_data')
end


ms_out = ms_in; % copy to maintain all other contents

% just for later check
segments_in = 1:length(ms_in.RawTraces);
keep_idx = segments_in(~ismember(segments_in, cells_to_remove));


%% remove all cells in the data fields with pre-segmented data
% fields_to_alter = {'time', 'RawTraces', 'FiltTraces', 'frameNum', 'vidNum'};

% for iF = fields_to_alter
for iC = sort(cells_to_remove, 'descend') % go backards or else everything is off and will not remove the right cells.
    ms_out.time(iC) = [];
    ms_out.RawTraces(iC) = [];
    ms_out.FiltTraces(iC) = [];
    ms_out.frameNum(iC) = [];
    ms_out.vidNum(iC) = [];
    
    fprintf('Removing cell: %0.f\n', iC);
end
% end



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
ms_out.format = [ms_out.format sprintf('; cell %.0f removed', cells_to_remove)];

if isfield(ms_out, 'history')
    ms_out.history{end+1} = {sprintf('; cell %.0f removed', cells_to_remove)};
else
    ms_out.history = {sprintf('; cell %.0f removed', cells_to_remove)};
end

