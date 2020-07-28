function events_out = MS_get_LFP_events_sandbox(cfg_in, csc)
%% MS_get_LFP_events_sandbox: WIP pipeline for rapid event detection in LFP data.
%
%
%
%    Inputs:
%     - cfg_in: [struct]    contains configuration (cfg)
%           see defaults in cfg_def below.  These will be overwritten by
%           the same field names in the input.
%
%     - csc [struct]        contains continously sampled channel (csc) data.  Should be
%           in the tsd format from MS_LoadCSC. Ideally this data should be
%           at 2000hz as this is what the built-in filters have been
%           optimized for.  Data with a higher sampling rate can be
%           decimated in the MS_LoadCSC function using
%           cfg.desired_sampling_frequency = 2000.
%
%
%
%    Outputs:
%     - events_out: [struct]   contains timestamps for the start and end of
%           each detected event as well as other user defined parameters
%           such as number of cycles, within event variance, min/max,...
%
%
%
%
% EC 2020-03-19   initial version based off of published methods from the
% van der Meer lab in Catanese et al. 2016, Carmichael et al., 2017,
% Carmichael et al. (2020; bioarxiv) and the vandermeer lab wiki:
% https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:dataanalysis
%
%%  Initialize

cfg_def = [];
cfg_def.check = 1; % plot all checks throughout.
%defaults to simple SWR detection.
cfg_def.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_def.filt.f  = [140 250]; % broad, could use 150-200?
cfg_def.filt.order = 4; %type filter order (fine for this f range)
cfg_def.filt.display_filter = 1; % use this to see the fvtool
cfg_def.filt.units = 'amplitude'; % units can also be in 'power'
cfg_def.nan_idx = [];


cfg = ProcessConfig(cfg_def, cfg_in);


%% basic filtering and thresholding
% mouse SWR parameters are based off of Liu, McAfee, & Heck 2017 https://www.nature.com/articles/s41598-017-09511-8#Sec6
%set up ripple band
%     cfg_filt_d = [];
%     cfg_filt_d.type = 'butter'; %Cheby1 is sharper than butter
%     cfg_filt_d.f  = [140 250]; % broad, could use 150-200?
%     cfg_filt_d.order = 4; %type filter order (fine for this f range)
%     cfg_filt_d.display_filter = 0; % use this to see the fvtool
% if cfg.check
%     cfg.filt.display_filter = 1; % use this to see the fvtool
% end
csc_filt = FilterLFP(cfg.filt, csc);


% convert to amplitude or power
amp_filt = csc_filt; % clone to make things simple and replace


for iChan = 1:size(csc_filt.data,1)
    amp_filt.data(iChan,:) = abs(hilbert(csc_filt.data(iChan,:)));
    % Convolve with a gaussian kernel (improves detection)
    if isfield(cfg, 'kernel')
        kernel = gausskernel(cfg.kernel.samples,cfg.kernel.sd); % note, units are in samples; for paper Methods, need to specify Gaussian SD in ms
        fprintf('\nGausskernal using %d samples = %0.0fms with SD = %d samples (%0.0fms)\n',cfg.kernel.samples, (cfg.kernel.samples/csc.cfg.hdr{1}.SamplingFrequency)*1000,cfg.kernel.sd,(cfg.kernel.sd/csc.cfg.hdr{1}.SamplingFrequency)*1000)
        amp_filt.data(iChan,:) = conv(amp_filt.data(iChan,:),kernel,'same');
    end
    amp_filt.units = 'amplitude';
    
    if isfield(cfg.filt, 'units') && strcmpi(cfg.filt.units, 'power')
        amp_filt.data(iChan,:) =  amp_filt.data(iChan,:).^2;
        amp_filt.units = 'power';
    end
end

if cfg.check
    figure(111)
    plot(csc.tvec, csc.data(1,:),'k',csc_filt.tvec, csc_filt.data(1,:), 'r',...
        amp_filt.tvec, amp_filt.data(1,:),'b')
    legend({'Raw', 'filt', 'Amp'})
end


%% remove large amplitude artifacts before SWR detection

if isfield(cfg, 'artif_det') && ~isempty(cfg.artif_det)
    csc_artif = csc;
    for iChan = 1:size(csc_filt.data,1)
        csc_artif.data(iChan,:) = abs(csc_artif.data(iChan,:)); % detect artifacts both ways
    end
    %
    %     switch cfg.artif_det.method
    %         case 'zscore'
    %         art_thresh = std(csc_artif.data(1,:))*cfg.artif_det.threshold;
    %     end
    
    
    artif_evts = TSDtoIV(cfg.artif_det,csc_artif);
    
    cfg_temp = [];
    cfg_temp.d = [-cfg.artif_det.rm_len cfg.artif_det.rm_len];
    artif_evts = ResizeIV(cfg_temp,artif_evts);
    
    
    % plot
    if cfg.check
        plot(113)
        cfg_plot=[];
        cfg_plot.display = 'iv'; % tsd, iv
        cfg_plot.target = csc.label{1};
        PlotTSDfromIV(cfg_plot,artif_evts,csc);
        %         hline(cfg.artif_det.threshold )
        pause(3); close all;
    end
    
    % zero pad artifacts to improve reliability of subsequent z-scoring
    artif_idx = TSD_getidx2(csc,artif_evts); % if error, try TSD_getidx (slower)
    for iChan = 1:size(csc_filt.data,1)
        csc_filt.data(iChan,artif_idx) = 0;
        amp_filt.data(iChan,artif_idx) = 0;
    end
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d large amplitude artifacts detected and zero-padded from csc_filt.\n',length(artif_evts.tstart));
    
end

% plot
if cfg.check
    %     plot(114)
    hold on
    plot(amp_filt.tvec, csc_filt.data(1,:),'g');
    %     hline(cfg.artif_det.threshold, '--g')
    plot(amp_filt.tvec, amp_filt.data(1,:),'-k');
    
    pause(3); close all;
end


%% isolate candidate events

% get the thresholds
cfg_detect = [];
cfg_detect.operation = '>';
cfg_detect.dcn = cfg_detect.operation; % b/c odd var naming in TSDtoIV
cfg_detect.method = cfg.method;
cfg_detect.threshold = cfg.threshold;
cfg_detect.target = csc.label{1};
cfg_detect.minlen = cfg.min_len; % mouse SWR: 40ms from Vandecasteele et al. 2015
cfg_detect.merge_thr = cfg.merge_thr; % merge events that are within 20ms of each other.
cfg_detect.bad_idx = cfg.nan_idx; % indices to remove before thresholding.
[evt_candidates,evt_thr] = MS_TSDtoIV(cfg_detect,amp_filt);

% % now apply to all data
% cfg_select = [];
% cfg_select.dcn = '>';
% cfg_select.method = 'raw';
% cfg_select.threshold = evt_thr;
% cfg_select.target = 'CSC1.ncs';
% cfg_select.minlen = cfg_detect.minlen;
%
% [evt_ids,~] = TSDtoIV(cfg_select,amp_filt);

%% compute general stats (cycle count, varience, filtered varience, mean, min, evt_len)
cfg_cc = [];
cfg_cc.threshold_type = 'raw';
cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
cfg_cc.filter_cfg = cfg.filt;
events_out = CountCycles(cfg_cc,csc,evt_candidates);

fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events detected initially.\n',length(evt_candidates.tstart));

if cfg.check
    cfg_plot = [];
    cfg_plot.display = 'iv';
    cfg_plot.mode = 'center';
    cfg_plot.width = 0.2;
    cfg_plot.target = csc.label{1};
    cfg_plot.title = 'var';
    PlotTSDfromIV(cfg_plot,events_out,csc);
    pause(2); close all;
end

%% exclude events with insufficient cycles - count how many exist above same threshold as used for detection


% get get the evetns with sufficient cycles.
if isfield(cfg, 'nCycles') && ~isempty(cfg.nCycles)
    cfg_Nc = [];
    if  isfield(cfg, 'nCycles_operation')
        cfg_Nc.operation = cfg.nCycles_operation;
    else
        cfg_Nc.operation = '>=';
    end
    cfg_Nc.threshold = cfg.nCycles;
    events_out = SelectIV(cfg_Nc,events_out,'nCycles');
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events remain after cycle count thresholding (%d cycle minimum).\n',length(events_out.tstart), cfg_Nc.threshold);
end
%% check for evnts that are too long.
% add in a user field for the length of the events (currently not used)
events_out.usr.evt_len = (events_out.tend - events_out.tstart)';

if isfield(cfg, 'max_len') && ~isempty(cfg.max_len)
    %     cfg_max_len = [];
    %     cfg_max_len.operation = '<';
    %     cfg_max_len.threshold = .1;
    events_out = SelectIV(cfg.max_len,events_out,'evt_len');
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>:: %d events remain after event length cutoff (%s %d ms removed).\n',length(events_out.tstart), cfg.max_len.operation, (cfg.max_len.threshold)*1000);
end


% if isfield(cfg, 'min_len') && ~isempty(cfg.min_len)
% %     cfg_max_len = [];
% %     cfg_max_len.operation = '<';
% %     cfg_max_len.threshold = .1;
%     events_out = SelectIV(cfg.min_len,events_out,'evt_len');
%
%     fprintf('\n<strong>MS_SWR_Ca2</strong>:: %d events remain after event length cutoff (%s %d ms removed).\n',length(events_out.tstart), cfg.min_len.operation, (cfg.min_len.threshold)*1000);
% end

%% check for evnts with high raw varience. 'var_raw' is added as a events_out.usr field in CountCycles
if isfield(cfg, 'var') && ~isempty(cfg.var)
    % cfg.var = [];
    % cfg.var.operation = '<';
    % cfg.var.threshold = 1;
    events_out = SelectIV(cfg.var,events_out,'var_raw');
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events remain after raw varience thresholding (''var_raw'' %s %.2f removed).\n',length(events_out.tstart),cfg.var.operation, cfg.var.threshold);
end


%% remove events that cooinside with artifacts.
if  isfield(cfg, 'artif_det') && ~isempty(cfg.artif_det)
    events_out = DifferenceIV([], events_out, artif_evts);
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events remain after removing those co-occuring with artifacts.\n',length(events_out.tstart));
end



%% check again
if cfg.check
    cfg_plot = [];
    cfg_plot.display = 'iv';
    cfg_plot.mode = 'center';
    cfg_plot.width = 0.2;
    cfg_plot.target = csc.label{1};
    cfg_plot.title = 'var_raw';
    PlotTSDfromIV(cfg_plot,events_out,csc);
    pause(3); close all;
end


