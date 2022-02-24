function ms_out = MS_restrict(ms,tstart, tend)
%% MS_restrict: restricts the input ms data structure between tstart and tend.
%
%
%
%    Inputs:
%    - ms: [struct] contains data in the 'ms' format with time and other
%    fields of the same length (ie: 'Binary', 'RawTraces', 'position',...)
%
%    - tstart: [nSample array]   time points to start
%
%    - tend: [nSample array]  time points to stop
%
%    Outputs:
%    -
%
%
%
%
% EC 2022-02-22  ported version
% based on restrict.m by:
% MvdM 2014-07-20 initial version
% youkitan 2016-11-27 edit: type checking update

if nargin <3
    error('Requires start/end times as either unique inputs ')
end


% check the format of the ms input.
if ~isfield(ms, 'time')
    error('ms input should be in the ms format with a ''time'' field')
end

ms_out = ms;

%% restrict

known_len = size(ms.time,1); % should always be the correct number of cells for the number of segments.

s_idx = nearest_idx3(tstart, ms.time);
e_idx = nearest_idx3(tend, ms.time);

fprintf('Resizing traces in between %.1fs to %.1fs (duration: %0.1fs)\n',ms.time(s_idx)/1000, ms.time(e_idx)/1000, (ms.time(e_idx) - ms.time(s_idx))/1000)

fields = fieldnames(ms);
for iF = 1:length(fields)
    if  ischar(ms.(fields{iF}))
        fprintf('Skipping in <strong>''%s''</strong>...\n', fields{iF})
        continue
    end
    
        
    
    field_size = size(ms.(fields{iF}));
    
    cell_idx = find(field_size == known_len,1);
    
    
    if ~isempty(cell_idx)
        fprintf('Resizing traces in <strong>''%s''</strong>...\n', fields{iF})
        
        if cell_idx ==1
            ms_out.(fields{iF}) = [];
            ms_out.(fields{iF}) = ms.(fields{iF})(s_idx:e_idx,:);
        elseif cell_idx ==2
            ms_out.(fields{iF}) = [];
            ms_out.(fields{iF}) = ms.(fields{iF})(:,s_idx:e_idx);
        end
    end
end
ms_out.numFrames = length(ms_out.time); 
ms_out.restricted = [tstart, tend];