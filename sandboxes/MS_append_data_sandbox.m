function ms_out = MS_append_data_sandbox(ms_in, varargin)
%% MS_append_data
%  MS_append_data will take in a miniscope data structure 'ms_in' and
%  addpend while performing some checks to ensure that everything lines up.
%
%   Example:
%       ms_out = MS_append_data(ms_in, 'csc1', csc_data)
%           This will add in a field to the ms structure with the name csc1
%           and the corresponding csc_data. If the data is in the from for
%           a cell, say for trial or presegmented data (mulitple recording
%           periods, or task/rest blocks, ...) it will check that the
%           number of cells is consistent with the number of segments
%           already present in ms_in for the time field.
%
%   Inputs:
%       - ms_in: [struct] the ms data structure from msGenerateVideoObj_ConcatenateBatch*
%       - variablename - value pairs: these should in the format
%           (ms_in,'variable_name1', variable value1,variable_name2', variable value2,... )
%
%   Outputs:
%       - ms_out: [struct] with appended data field
%
%
%  * in https://github.com/torilynntemple/MiniscopeAnalysisConcatenation by
%   Tori-Lynn Temple (Williams lab at McGill)
%
%   EC 2019-12-13
%       - initial sandbox with basic functions
%% initialize

% get additional variable - value pairs and hold them in the workplace by
% their variable name.
if nargin
    for iV = 1:2:length(varargin)
        eval([varargin{iV}, ' = ', 'varargin{iV+1};']);
        var_names{iV} = varargin{iV}; % list of variable names
        var_val{iV} = varargin{iV+1}; % variable values
    end
end

% clone the ms_in to maintain it for output.
ms_out = ms_in;


%% append the data w/ checks

for iV = 1:length(var_names)
    
    if ~iscell(var_val{iV})
        ms_out.(var_names{iV}) = struct(var_names{iV}, var_val{iV});
    else
        % run checks on the number of segments in the ms_in.times vs the
        % input variable
        if iscell(ms_in.time) && length(var_val{iV}) ~= length(ms_in.time)
            error('input variable ''%s'' is not the same size as ms_in.time and should be examined to avoid misaligned data', var_names{iV})
        else
            ms_out.(var_names{iV}) = struct(var_names{iV}, var_val{iV}); % if they are the same size than appended them.
        end
    end
end

%% clean up

if isfield(ms_out, 'history')
    for ii = 1:length(var_names)
        ms_out.history{end+1} = {sprintf('MS_append_data on %s  added ''%s''', date, var_names{ii})};
    end
else
    for ii = 1:length(var_names)
        ms_out.history = {sprintf('MS_append_data on %s  added ''%s''', date, var_names{ii})};
    end
end

