function [data_out] = HF_preprocess(phy2_dir, csc_dir, evts_dir)
% HF_preprocess loads and preprocesses all the principle data types for a
% recording session. 

%% load the phy along with the params file
params = OE_load_params(phy2_dir);

S = OE_phy2TS(phy2_dir);


evts = OE_load_binary_evts(evts_dir);
%% load the lfp
csc_list = dir([csc_dir filesep '*CH*.continuous']);

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
    adc_ts.t{1} = adc_ts.t{1} - adc_f_tsd.tvec(1); 
    rate_tsd.tvec = rate_tsd.tvec - rate_tsd.tvec(1); 
end
%% correct the events times and the LFP so that it aligns with the spikes

% correct for start of csc. Why is this offset? 
csc.tvec = csc.tvec- csc.tvec(1); 
% loop over events and remove the offset. 
for ii = 1:length(evts.t)
    evts.t{ii} = evts.t{ii} - offset; 
end

% compile the data as one object

data_out.params = params; 
data_out.S = S;
data_out.csc = csc; 
data_out.evts = evts; 


