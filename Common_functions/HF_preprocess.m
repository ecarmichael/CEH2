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

S = OE_phy2TS(phy2_dir);


evts = OE_load_binary_evts(evts_dir);
%% load the lfp
csc_list = dir([csc_dir filesep '*CH*.continuous']);

if ~isempty(csc_idx)
    csc_list(~(1:length(csc_list) == csc_idx))= []; 
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

csc.data = csc.data'; 
csc.label = labels; 

offset = csc.tvec(1); 

cfg_in.decimateFactor = 15; 
csc = decimate_tsd(cfg_in, csc);


%% get the AD channel piezo
adc = dir([csc_dir filesep '*ADC*.continuous']); 
if ~isempty(adc)
    [adc_ts, adc_f_tsd, rate_tsd] =  HF_piezo2ts(adc(1).name, 0); 
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


% compile the data as one object

data_out.params = params; 
data_out.S = S;
data_out.csc = csc; 
data_out.evts = evts; 
data_out.licks = rate_tsd; 
if ~isempty(vr_fname)
    data_out.vr = vr; 
end


end


