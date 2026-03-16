function csc = MS_LoadCSC_CE(cfg_in)



cfg_def = [];
cfg_def.fc = {};
cfg_def.decimateByFactor = [];
cfg_def.desired_sampling_frequency = 2000;
cfg = ProcessConfig(cfg_def, cfg_in);




csc_list = dir(['*CH*.continuous']);

% sort the csc based on channel number. 
for ii = length(csc_list):-1:1
    ch_idx =  strfind(csc_list(ii).name,'_CH'); 
    con_idx =  strfind(csc_list(ii).name,'.continuous'); 
    csc_num(ii) = str2double(csc_list(ii).name(ch_idx+3:con_idx)); 
    csc_name{ii} = csc_list(ii).name; 
end

% sort
[~, sort_idx] = sort(csc_num); 
csc_list = csc_list(sort_idx); 

if  ~isempty(cfg.fc)  && iscell(cfg.fc) % if csc_idx is a string look for those patterns
    
    temp_list = [];

    for ii = 1:length(cfg.fc)
        temp_idx = find(contains(csc_name, [cfg.fc{ii} '.continuous']));
        if length(temp_idx) > 1
            error('Too many csc files containing this name')
        end
        if ~isempty(temp_idx)
            temp_list(ii).name= csc_name{temp_idx};
            temp_list(ii).folder= csc_list(1).folder;
        end
    end

    csc_list = temp_list;

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

csc.cfg.hdr{1}.og_sampleRate = csc.cfg.hdr{1}.sampleRate; 
csc.cfg.hdr{1} = rmfield(csc.cfg.hdr{1}, 'sampleRate'); 

% csc.tvec = ;
% csc.tvec = csc.tvec - csc.tvec(1) + (csc.tvec(1) - ts_prime(1)); % zero out the csc. 
if ~isempty(cfg.desired_sampling_frequency)

        if mod(csc.cfg.hdr{1}.SamplingFrequency, cfg.desired_sampling_frequency) ~=0
            error('Current Sampling Frequency cannot be divided by desired Sampling Frequency')
        end
        
        this_Fs = csc.cfg.hdr{1}.SamplingFrequency;
        this_decimate = this_Fs / cfg.desired_sampling_frequency; % must be a whole number.
        fprintf('%s: Current Sampling Frequency: %dHz...\n',mfilename,this_Fs);


cfg_decimate = [];
cfg_decimate.decimateFactor = this_decimate; 
csc = decimate_tsd(cfg_decimate, csc);
end