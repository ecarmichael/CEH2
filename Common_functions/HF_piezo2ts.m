function [adc_ts, adc_f_tsd, rate_tsd] =  HF_piezo2ts(fname, plt_flg)
%% HF_piezo2iv: loads the Analog signal from the piezo, applies some filters and extracts the deflections as timestamps.
%
%  to convert to a rate use MS_spike2rate.m
%
%

% Initialize 
if nargin   < 2
    plt_flg = true; % Set default plotting flag if not provided
end


%% load the data into the tsd format.

[data, tvec, info] = load_open_ephys_data(fname);
adc.type = 'tsd';
adc.units = 'au';
adc.tvec = tvec;
adc.data = data';
adc.label = {'Piezo'};
adc.cfg.hdr{1} = info.header;
adc.cfg.hdr{1}.SamplingFrequency = info.header.sampleRate;
adc.cfg.history.mfun{1} = 'tsd';
adc.cfg.history.cfg{1} = [];


adc = MS_decimate_CSC(adc, adc.cfg.hdr{1}.SamplingFrequency./500);
%% filter the signal 


cfg_filt = [];
cfg_filt.f = [1 55];  % seems to capture both slow and fast well
cfg_filt.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt.order = 6; %type filter order
adc_f_tsd = FilterLFP(cfg_filt, adc);


%% find peaks and convert to licks
[~, l_idx] = findpeaks(abs(adc_f_tsd.data), 'MinPeakHeight',0.005, 'MinPeakDistance', adc_f_tsd.cfg.hdr{1}.SamplingFrequency*0.1); 

fprintf('Lick rate = %0.2f/s\n', length(l_idx)./(adc.tvec(end) - adc.tvec(1)))

if plt_flg
    figure
    clf
    findpeaks(abs(adc_f_tsd.data), 'MinPeakHeight',0.01, 'MinPeakDistance', adc_f_tsd.cfg.hdr{1}.SamplingFrequency*0.1)
    hold on
    plot(adc.data, 'color', 'black')
end


% adc_ts = ts(adc.tvec(l_idx), {'licks'})

adc_ts.type = 'ts';
adc_ts.t = {adc.tvec(l_idx)};
adc_ts.label = {'licks'};

rate_tsd = MS_spike2rate(adc_ts, adc_f_tsd.tvec, .05, .1);


