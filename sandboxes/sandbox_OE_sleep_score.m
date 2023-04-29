[data, tvec, hdr] = load_open_ephys_data('100_RhythmData-B_CH9.continuous');

cfg.decimateByFactor = 30000 / 1000;
data = decimate(data,cfg.decimateByFactor);
tvec = tvec(1:cfg.decimateByFactor:end);

csc = tsd(tvec, data', 'CSC7.ncs');
csc.cfg.hdr{1}.SamplingFrequency = 1000;
clear data tvec
%%
[emg_data, emg_tvec, emg_hdr] = load_open_ephys_data('100_RhythmData-B_CH52.continuous');

emg_data = decimate(emg_data,cfg.decimateByFactor);
emg_tvec = emg_tvec(1:cfg.decimateByFactor:end);

emg = tsd(emg_tvec, emg_data', 'EMG');
emg.cfg.hdr{1}.SamplingFrequency = 1000;
clear emg_data emg_tvec

%% get the theta delta ratio
cfg_emg = [];
cfg_emg.threshold = 0;
cfg_emg.f = [15 30]; 
% cfg_emg.type = 'cheby1';
% cfg_emg.display_filter = 1
emg_f = FilterLFP(cfg_emg,csc);
%% sleep state
wake_t = [0 1660 11177 13440 17000 csc.tvec(end)]; 

wake_idx = nearest_idx(wake_t, csc.tvec);
wake_idx = reshape(wake_idx,2, length(wake_idx)/2)'; 
[hypno, csc, emg] = dSub_Sleep_screener(csc, emg_f, wake_idx)