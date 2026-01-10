function [SWR_evts, cfg_swr] = MS_SWR_detector(csc, chan_to_use, plot_flag)
%% MS_SWR_detector: default SWR detector. Meant for speed, not accuracy.
%
%
%
%    Inputs:
%    - csc: [struct]      data in the TSD format. Should one have one data
%    channel.
%
%    - chan_to_use: [sting]  Label for the channel to use.
%
%    Outputs:
%    - SWR_evts: [struct] SWR times in the TS format.
%
%
%
%
% EC 2023-06-06   initial version
%
%
%
%% initialize

if nargin == 1
    chan_to_use = [];
    plot_flag = 1;
elseif nargin == 2
    plot_flag = 1;
    
end


if length(csc.label) >1 && isempty(chan_to_use);
    error('Too many data channels. Filtering and detection is only designed for one channel ATM')
end

if ~isempty(chan_to_use)
    chan_idx = find(strcmp(csc.label, chan_to_use));
    
    temp_data = csc.data(chan_idx,:);
    csc.data = [];
    csc.data(1,:) = temp_data;
    temp_label = csc.label{chan_idx};
    csc.label = [];
    csc.label{1} = temp_label;
    
    temp_hdr = csc.cfg.hdr{chan_idx};
    csc.cfg.hdr = [];
    csc.cfg.hdr{1} = temp_hdr;
    
    clear temp_*
end

%% filter the LFP into the ripple band.
cfg_swr = [];
cfg_swr.check = 0; % plot checks.
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [125 250]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter = 0; % use this to see the fvtool

% smoothing
cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;

% detection
cfg_swr.artif_det.method = 'zscore';
cfg_swr.artif_det.threshold = 5;
cfg_swr.artif_det.dcn = '>';
cfg_swr.artif_det.rm_len = .2;
cfg_swr.threshold = 1.5;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
cfg_swr.merge_thr = 0.02; %merge events that are within 20ms of each other.

%
% cfg_swr.nan_idx = nan_idx; % where are any nans, say from excluding artifacts, other events...

% restrictions
cfg_swr.max_len = [];
cfg_swr.max_len.operation = '<';
cfg_swr.max_len.threshold = .1;

%                 cfg_swr.min_len = [];
%                 cfg_swr.min_len.operation = '<';
%                 cfg_swr.min_len.threshold = .2;
cfg_swr.nCycles = 5; % number of cycles
cfg_swr.nCycles_operation = '>='; % number of cycles

% variaence
% cfg_swr.var = [];
% cfg_swr.var.operation = '<';
% cfg_swr.var.threshold = 1;

[SWR_evts, csc_filt] = MS_get_LFP_events_sandbox(cfg_swr, csc);

%
%         cfg_max_len = [];
%         cfg_max_len.operation = '>';
%         cfg_max_len.threshold = 5;
%          SWR_evts = SelectIV(cfg_max_len,SWR_evts,'nCycles');
% check quality.
if plot_flag
    cfg_plot.display = 'iv'; %'iv';
    cfg_plot.title = 'var';
    PlotTSDfromIV(cfg_plot, SWR_evts, csc)
    
end

%% figure(909) SWR triggered average
if plot_flag
    swr_cent = IVcenters(SWR_evts);
    win = csc.cfg.hdr{1}.SamplingFrequency / 2;
    
    swr_idx = nearest_idx3(swr_cent, csc.tvec);
    
    all_SWR = [];
    for ii = length(swr_cent):-1:1
        if (swr_idx(ii) + win) < length(csc.tvec)
            all_SWR(ii,:) = csc.data(1,swr_idx(ii) - win: swr_idx(ii)+win);
        end
    end
    
    
    figure(909)
    plot(-.5:1/(csc.cfg.hdr{1}.SamplingFrequency):.5, mean(all_SWR))
    
end
