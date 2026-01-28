function [data_out] = HF_preprocess(phy2_dir, csc_dir, evts_dir, vr_fname, csc_idx)
% HF_preprocess loads and preprocesses all the principle data types for a
% recording session.

if nargin < 4
    csc_idx = [];
    vr_fname  =[];
elseif nargin < 5
    csc_idx = [];
end
%% load the phy along with the params file
params = OE_load_params(phy2_dir);

if ~isempty(phy2_dir)
S = OE_phy2TS(phy2_dir, params);
else
    S = []; 
end


evts = OE_load_binary_evts(evts_dir);

OE_rec = readNPY([phy2_dir filesep 'timestamps.npy']);
OE_rec = [OE_rec(1) OE_rec(end)];

% align to the recording start time. (unclear why this is not the same as
% the spike times. 
S_r = restrict(S, rec_iv);

for ii = 1:length(S_r.t)
    S_r.t{ii} = S_r.t{ii} - rec_iv.tstart;
end
%% load the lfp
csc_list = dir([csc_dir filesep '*CH*.continuous']);

% sort the csc based on channel number. 
for ii = length(csc_list):-1:1
    ch_idx =  strfind(csc_list(ii).name,'_CH'); 
    con_idx =  strfind(csc_list(ii).name,'.continuous'); 
    csc_num(ii) = str2double(csc_list(ii).name(ch_idx+3:con_idx)); 
end

% sort
[~, sort_idx] = sort(csc_num); 
csc_list = csc_list(sort_idx); 


if ~isempty(csc_idx)
    csc_list(~ismember(1:length(csc_list), csc_idx))= [];
end

csc= []; labels = [];
for ii = 1:length(csc_list)

    if ii == 1
        [data, tvec, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc = tsd(tvec, data);
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;

    else
        [data, ~, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc.data =[csc.data, data];
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;
    end
end
fs = csc.cfg.hdr{1}.SamplingFrequency; 

csc.data = csc.data';
csc.label = labels;

offset = csc.tvec(1);

cfg_in.decimateFactor = 15;
csc = decimate_tsd(cfg_in, csc);

% load the OE version of the events.
evts_list = dir([csc_dir filesep '*Data*.events']);

OE_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);
%% get the AD channel piezo
adc = dir([csc_dir filesep '*ADC*.continuous']);
if ~isempty(adc)
    [adc_ts, adc_f_tsd, rate_tsd] =  HF_piezo2ts([adc(1).folder filesep adc(1).name], 0);
    adc_f_tsd.tvec = adc_f_tsd.tvec - adc_f_tsd.tvec(1);
    rate_tsd.tvec = rate_tsd.tvec - rate_tsd.tvec(1);

    evts.t{end+1} = adc_ts.t{1};
    evts.label{end+1} = 'Licks';
else
    rate_tsd = tsd();
end


%% correct the events times and the LFP so that it aligns with the spikes

% correct for start of csc. Why is this offset?
csc.tvec = csc.tvec- csc.tvec(1);
% loop over events and remove the offset.
for ii = 1:length(evts.t)
    evts.t{ii} = evts.t{ii} - offset;
end

% loop over spikes and remove the offset.
for ii = 1:length(evts.t)
    evts.t{ii} = evts.t{ii} - offset;
end


% loop over events and remove the offset.
for ii = 1:length(OE_evts.t)
    OE_evts.t{ii} = OE_evts.t{ii} - offset;
end


%% load the log file
if ~isempty(vr_fname)
    vr = HF_load_VR(vr_fname);

    % align to the first reward event
    first_reward_oe = evts.t{ismember(evts.label, '4')}(1,1);

    % first vr reward
    first_reward_vr = vr.evt.t{contains(vr.evt.label, 'Collision with Rwd1')}(1);
    vr_start = first_reward_vr - vr.pos.tvec(1);

    % Align the OE data to the first reward event
    vr.pos.tvec = vr.pos.tvec - vr_start + first_reward_oe;

    for ii = 1:length(vr.evt.t)

        vr.evt.t{ii} = vr.evt.t{ii} - vr_start + first_reward_oe;


    end

end

% compile the data as one object

data_out.params = params;
data_out.S = S;
data_out.csc = csc;
data_out.evts = evts;
data_out.OE_evts = OE_evts;
data_out.licks = rate_tsd;
data_out.licks.data(end+1,:) = adc_f_tsd.data(1:end-1); 
data_out.licks.label{end+1} = 'Filt ADC'; 
if ~isempty(vr_fname)
    data_out.vr = vr;
end





