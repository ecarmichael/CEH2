function ds_iv = HF_dentate_detect(data_in)
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

csc_f = FilterLFP(cfg_filter, this_csc);

d_fs = floor(1/(50/csc_f.cfg.hdr{1}.SamplingFrequency)); % 50ms window. 


%% detect



%% classify


