function csc = OE_old_csc2TSD(cfg_in)
%% OE_old_csc2mat.m: loads and converts OE file types into TSD
%
%
%
%    Inputs:
%    -
%
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2023-06-05   initial version
%
%
%
%% initialize


cfg_def = [];
cfg_def.fc = {};
cfg_def.decimateByFactor = [];
cfg_def.desired_sampling_frequency = 2000;
cfg = ProcessConfig(cfg_def, cfg_in);


%% get all the files in a dir

if isempty(cfg.fc)
    fnames = dir('*.continuous');
    for ii = 1:length(fnames)
        
        cfg.fc{ii} = fnames(ii).name;
        
        [parts]= strsplit(cfg.fc{ii}, '.');
        [label_num{ii}]= parts{1}(strfind(parts{1}, 'CH'):end);
        
        this_num = regexp(label_num{ii},'\d*','Match');
        num_labels(ii) = str2num(this_num{1});
        
        
    end
    
    [~, sort_idx] = sort(num_labels);
    cfg.fc =  cfg.fc(sort_idx);
    
end


%% Load the OE data and save it for TSD
data = [];
tvec = [];
hdr = [];

for ii = length(cfg.fc):-1:1
    
    [this_data, this_tvec, this_hdr] = load_open_ephys_data(cfg.fc{ii});
    
    % decimate data if specified
    if ~isempty(cfg.decimateByFactor) && isempty(cfg.desired_sampling_frequency) % make sure you are not trying to decimate and not also automating the downsampling.
        
        fprintf('%s: Decimating by factor %d...\n',mfun,cfg.decimateByFactor)
        this_data = decimate(this_data,cfg.decimateByFactor);
        this_tvec = this_tvec(1:cfg.decimateByFactor:end);
        this_hdr.SamplingFrequency = this_hdr.header.sampleRate./cfg.decimateByFactor;
        
    end
    
    % "automated decimation" (get a desired sampling frequency given the
    % data). Save the time/memory from having to do this after.
    % decimate data if specified
    if ~isempty(cfg.desired_sampling_frequency) && isempty(cfg.decimateByFactor) && (this_hdr.header.sampleRate~= cfg.desired_sampling_frequency) % make sure you are not trying to decimate and not also automating the downsampling.
        
        if mod(this_hdr.header.sampleRate, cfg.desired_sampling_frequency) ~=0
            error('Current Sampling Frequency cannot be divided by desired Sampling Frequency')
        end
        
        this_Fs = this_hdr.header.sampleRate;
        this_decimate = this_Fs / cfg.desired_sampling_frequency; % must be a whole number.
        fprintf('%s: Current Sampling Frequency: %dHz...\n',mfilename,this_Fs);
        
        
        this_data = decimate(this_data,this_decimate);
        this_tvec = this_tvec(1:this_decimate:end);
        this_hdr.SamplingFrequency = this_hdr.header.sampleRate./this_decimate;
        this_hdr.header.sampleRate = this_hdr.header.sampleRate./this_decimate;
        fprintf('%s: Decimating by factor %d to output data at desired Fs of %dHz \n',mfilename,this_decimate, cfg.desired_sampling_frequency)
        
        
    end
    
    
    [parts]= strsplit(cfg.fc{ii}, '.');
    
    [labels{ii}]= parts{1}(strfind(parts{1}, 'CH'):end);
    
    data(ii,:) = this_data;
    tvec(ii, :) = this_tvec;
    hdr{ii} = this_hdr;
end

%% convert to TSD


csc = tsd(tvec(1,:)', data, labels);

csc.cfg.hdr = hdr;
% csc.cfg.hdr.SamplingFrequency = hdr{1}.header.sampleRate;

for ii = 1:length(csc.label)
    % csc.cfg.hdr{ii}.SamplingFrequency = round(1/mode(diff(csc.tvec)));
    csc.cfg.hdr{ii}.SamplingFrequency = hdr{1}.header.sampleRate;
    
end







