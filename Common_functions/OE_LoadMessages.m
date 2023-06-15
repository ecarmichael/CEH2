function evt_out = OE_LoadMessages(fname)
%% OE_LoadMessages: load the messages.events file from OE
%
%
%
%    Inputs:
%    - fname: [string]      name of the file to load. If empty it looks for
%    messages.events in the current dir.

%
%
%    Outputs:
%    - evts: [struct]      event times and messages in the TS format.
%
%
%
%
% EC 2023-06-08   initial version
%
%
%
%% initialize
if nargin < 1
    fname  = 'messages.events';
end 

%filetype = fname(max(strfind(fname,'.'))+1:end); % parse filetype

fid = fopen(fname);

fseek(fid,0,'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');


NUM_HEADER_BYTES = 1024;


disp(['Loading events file...']);

index = 0;

hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
hdr_char = char(hdr')

n_lines = regexp(hdr_char, '[\n]');

ts_in = []; label = [];
for ii =length(n_lines):-1:1
    
    if ii == 1
        parts = strsplit(hdr_char(1:n_lines(ii)-1), ',');
        
    else
        parts = strsplit(hdr_char(n_lines(ii-1):n_lines(ii)-1), ',');
        
    end
    
    ts_ret = regexp(parts{1}, '[\n]') ;
    
    if ~isempty(ts_ret)
        ts_in{ii} = str2double(parts{1}(ts_ret+1:end));
    else
        ts_in{ii} = str2double(parts{1});
    end
    label{ii} = parts{end}(2:end);
    
end

% %% extract the sampling frequency and normalize the values. 

for ii = 1:length(label)
    
    at_idx = strfind(label{ii}, '@');
    hz_idx = strfind(label{ii}, 'Hz');
    
    if ~isempty(at_idx) && ~isempty(hz_idx)
        this_num = regexp(label{ii}(at_idx:hz_idx),'\d*','Match');
        Fs= str2double(this_num{1});
    end
end

if length(Fs) ~=1
    error('Could not detect Sampling frequency in the messages')
end


for ii = length(ts_in):-1:2
   ts_in{ii} = ts_in{ii}./Fs;  
end

evt_out = ts(ts_in);
evt_out.label = label;
evt_out.cfg.hdr{1}.SamplingFrequency = Fs; 

