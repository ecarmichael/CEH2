function ds_iv = HF_dentate_detect(data_in, dg_chan, swr_iv)
%% HF_dentate_detect:
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
% EC 2026-11-10   initial version 
%
%
%
%% initialize
if nargin < 3
    swr_iv =[]; 
end


%% filter:  LDS et al. 2026
% sleep/rest periods
% LFPs from the dentate gyrus region were initially referenced by subtracting 
% signals from a channel located at the top of the probe [done in recording with virtual reference]
% filtered within a frequency range of 1–200 Hz/ 4th-order Butterworth filter.
% Peaks seven times the median absolute value were identified as candidate dentate spike events

cfg_filter.type = 'butter'; %Cheby1 is sharper than butter
cfg_filter.f  = [1 200]; % broad, could use 150-200?
cfg_filter.order = 4; %type filter order (fine for this f range)
cfg_filter.display_filter = 0; % use this to see the fvtool

csc_f = FilterLFP(cfg_filter, data_in);

d_fs = floor(1/(50/csc_f.cfg.hdr{1}.SamplingFrequency)); % 50ms window. 

amp_filt = abs(hilbert(csc_f.data)); 
%% quick plot to check filter

figure(1010)
clf

hold on
for ii = 1:size(csc_f.data,1)

    if ii == dg_chan
        plot(data_in.tvec, (data_in.data(ii,:)/1000)+ii, 'LineWidth', 2)
        plot(data_in.tvec, (csc_f.data(ii,:)/1000)+ii, 'k','LineWidth', 2)
        plot(data_in.tvec, (amp_filt(ii,:)/1000)+ii, 'k','LineWidth', 2)

    else
        plot(data_in.tvec, (data_in.data(ii,:)/1000)+ii)
        plot(data_in.tvec, (csc_f.data(ii,:)/1000)+ii, 'k')
        plot(data_in.tvec, (amp_filt(ii,:)/1000)+ii, 'k')
    end
end

xlim([200 300])

set(gca, 'YDir', 'reverse')

if ~isempty(swr_iv)
    xline(IVcenters(swr_iv)); 
end


%% detect

cfg_detect = [];
cfg_detect.operation = '>';
cfg_detect.dcn = cfg_detect.operation; % b/c odd var naming in TSDtoIV
cfg_detect.method = 'raw';
cfg_detect.threshold = 7;
cfg_detect.target = data_in.label{dg_chan};
cfg_detect.minlen = 0.05; % mouse SWR: 40ms from Vandecasteele et al. 2015
cfg_detect.merge_thr = 30; % merge events that are within 20ms of each other.
[evt_candidates,evt_thr] = MS_TSDtoIV(cfg_detect,amp_filt);



%% classify


