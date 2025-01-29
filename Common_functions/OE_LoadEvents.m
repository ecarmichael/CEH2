function Evt_out = OE_LoadEvents(fname, Fs)
%% OE_LoadEvents: load the *.events file from OE
%
%
%
%    Inputs:
%    - fname: [string]      name of the file to load. If empty it looks for
%    messages.events in the current dir.
%
%   - Fs: [double]         sampling frequency in the system. Needed to
%    convert timestamps into time. Default assumes 30,000 hz. 
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

if nargin <1
    
    files  = dir('*.events');
    fprintf('<strong>%s</strong>: %d events file detected...', mfilename, length(files))

    for ii = 1:length(files)
        if ~contains( files(ii).name, 'messages')
            fname = files(ii).name;
            fprintf('using: <strong>%s</strong>\n', fname)
        end
    end
    
    Fs = 30000; 
    disp('No Fs specified. Assuming 30kHz...')
elseif nargin ==2
    Fs = 30000; 
        disp('No Fs specified. Assuming 30kHz...')
end

%% load the events file

[evts_label, evts_ts, evts_hdr] = load_open_ephys_data(fname);

evts_ts = evts_ts ./ Fs; % correct samples to time. 

%% get all of the unique entries

u_labels = unique(evts_label); 
t = [];
label = [];
for ii = length(u_labels):-1:1
   
    this_idx =  evts_label ==  u_labels(ii);
    
    t{ii} = evts_ts(this_idx);
    label{ii} = num2str(u_labels(ii));
    
end

% convert to TS format

Evt_out = [];

Evt_out.type= 'ts';
Evt_out.t = t;
Evt_out.label= label;
Evt_out.cfg.history.mfun{1} = mfilename;
Evt_out.cfg.history.cfg{1} = [];
Evt_out.cfg.hdr{1} = evts_hdr; 


