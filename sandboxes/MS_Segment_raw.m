function ms_seg_resize = MS_Segment_raw(cfg_in, csc_dir, ms_dir, ms_resize_dir)
%% MS_Segment_raw: Load, segment, visualize, and resize miniscope ('ms') and neuralynx ('nlx') data.
%
%
%
%    Inputs:
%     - cfg_in: [struct] contains user configurations that wil overwrite
%     the defaults in cfg_def.
%
%     - csc_dir: [string] path to nlx data (contains csc.ncs and Events.nev
%     files)
%
%     - ms_dir: [string] path to the miniscope data folders containing the
%     ms.mat file and subfolders for each recording block.
%
%     - ms_resize_dir: [string] [optional] path to where the resized
%     ms_resize.mat file should be saved.  If this input is empty it will
%     be the same as ms_dir.
%
%    Outputs:
%     - ms_resize: [struct] the segmented and resized ms structure with the
%     segmented nlx csc and evt data structures appended to it. toggle
%     cfg.append_nlx to remove the nlx data from ms_resize before saving.
%
%
%
%
% EC 2020-02-18   initial version
%
%
%% initialize

global PARAMS
if nargin < 2
    error('Insufficient inputs.  Requires at least: cfg_in, csc_dir')
elseif nargin ==2
    ms_dir = csc_dir;
    ms_resize_dir = ms_dir;
    warning('No ''ms_dir'' OR ''ms_resize_dir'' specified, using csc_dir for both...')
elseif nargin ==3
    ms_resize_dir = ms_dir;
    warning('No ''ms_resize_dir'' specified, using ms_dir...')
end

fprintf('\n<strong>MS_Segment_raw</strong>: csc_dir = %s \n', csc_dir);
fprintf('\n                ms_dir = %s \n', ms_dir);
fprintf('\n                ms_resize_dir = %s \n', ms_resize_dir);


% set the default configs.
cfg_def = [];
cfg_def.append_nlx = 1; % can be 0 to remove the nlx data structs from ms_resize before saving/outputting
cfg_def.save_ms_resize = 1; % can be 0 to not save the ms_resize and instead just take the output.
cfg_def.csc.fc = {'CSC1.ncs','CSC7.ncs'}; % use csc files from Keys. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
cfg_def.csc.label = {'EMG', 'LFP'}; % custom naming for each channel.
cfg_def.csc.desired_sampling_frequency = 2000;
cfg_def.evt.comb_chans = [3 4]; % see below.
cfg_def.evt.t_chan = 5; % which Events channel to use.  Seems to be 3 (ms TTL on) 4(ms TTL off) which combined make a new channel 5

% filters
cfg_def.filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_def.filt_d.f  = [1 5];
cfg_def.filt_d.order = 8; %type filter order
cfg_def.filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_def.filt_t.f  = [6 11];
cfg_def.filt_t.order = 3; %type filter order

% matching and resizing
cfg_def.TS_nlx_match = 1; % use this if the TS file is off by one sample. It will trim the ms data to fit the nlx data. [only a bandaid solution until segmentation can be verified with ground truth]
cfg_def.resize = 1; % use this to resize the data using the selection tool.
cfg_def.spec.win_s = 2^10; % window size for spec keep in base 2


cfg = ProcessConfig(cfg_def, cfg_in);


%% load the Miniscope data
% move to the ms directory
cd(ms_dir)

%load the ms.mat file.
load('ms.mat')


[TS, TS_name] = MS_collect_timestamps(ms_dir);

% get the hypnogram labels
[hypno_labels, time_labels] = MS_get_hypno_label([], TS_name);

% compare to TS to ms
fprintf('\n****Comparing TS files to processed miniscope (ms) data\n')
for iT = 1:length(TS)
    if length(TS{iT}.system_clock{1}) == ms.timestamps(iT)
        disp([TS_name{iT}   ':  ' num2str(length(TS{iT}.system_clock{1}))   ' - ms TS: ' num2str(ms.timestamps(iT))  '   ~ ' num2str(length(TS{iT}.system_clock{1}) / TS{iT}.cfg.Fs{1}) 's'])
    else
        warning(['TS do not match ms data' TS{iT}.filename   ':  ' num2str(length(TS{iT}.system_clock{1}))   ' - ms TS: ' num2str(ms.timestamps(iT))])
    end
end
fprintf('\n<strong>MS_Segment_raw</strong>: ms.mat loaded and timestamps collected\n');

%% Get some basic statistics for the individual Ca Traces
cfg_stats = [];
cfg_stats.data_type = 'RawTraces';

ms = MS_characterize_trace(cfg_stats, ms);
fprintf('\n<strong>MS_Segment_raw</strong>: basic statstics computed for each Ca trace\n');


%% append the TS folder name to the ms struct.

ms = MS_append_data_sandbox(ms, 'file_names', TS_name);


%% remove cells with low quality signals. [ can be used to remove cells with specific statistics in the raw traces.
% cfg_remove_trace = [];
% cfg_remove_trace.threshold = 8;
% cfg_remove_trace.remove_idx = ms.stats.RawTraces.z_max <= cfg_remove_trace.threshold;
%
% ms = MS_Remove_trace(cfg_remove_trace, ms);
% fprintf('\n<strong>MS_Segment_raw</strong>: %d Traces had max zscore values  < %d. Taken to indicate poor spike quality \n',sum(cfg_remove_trace.remove_idx), cfg_remove_trace.threshold);



%%  convert the Ca transitents into a binarized vector. Can be used to binarize the data using peaks and smoothing.
% cfg_bin = [];
% cfg_bin.method = 'rise';
% cfg_bin.rise.smooth_type = 'sgolay';
% cfg_bin.threshold = 2;
% ms = MS_binarize_data_sandbox(cfg_bin, ms);
% fprintf('\n<strong>MS_Segment_raw</strong>: miniscope data has been binarized using a %s method with a threshold of %d\n', cfg_bin.method, cfg_bin.threshold);
% pause(2)
% close all


%% segment the data
cfg_seg = [];
% cfg_seg.user_fields = {'BinaryTraces'};
ms_seg = MS_segment_ms_sandbox(cfg_seg, ms);

fprintf('\n<strong>MS_Segment_raw</strong>: miniscope data has been segmented into %d individual recording epochs\n method used: %s\n', length(ms_seg.time), ms_seg.format);

%% Load nlx data
cd(csc_dir)
% load the Keys file with all of the experiment details.
%(can be generated with the 'MS_Write_Keys' function)
if exist('*Keys.m', 'file')
    ExpKeys = MS_Load_Keys();
end

% load events
nlx_evts = LoadEvents([]); % get '.nev' file in this dir.

% load the NLX CSC data (using vandermeer lab code) [todo:replace with own]

csc = MS_LoadCSC(cfg.csc); % need to comment out ExpKeys lines in LoadCSC


% extract NLX event epochs

nlx_evts.t{cfg.evt.t_chan} = sort([nlx_evts.t{cfg.evt.comb_chans(1)}, nlx_evts.t{cfg.evt.comb_chans(2)}]);
nlx_evts.label{5} = ['merge TTls at ' num2str(cfg.evt.comb_chans(1)) ' and ' num2str(cfg.evt.comb_chans(2))];
[evt_blocks, ~, evt_duration] = MS_extract_NLX_blocks_sandbox(cfg.evt, nlx_evts);
pause(1)
close

% compare to TS to ms
fprintf('\n****Comparing TS files to processed miniscope (ms) data\n')
if length(TS) ~= length(evt_blocks)
    warning('Number of Timestamp files (%s) does not match the number of detected NLX event blocks (%s)', num2str(length(TS)),num2str(length(evt_blocks)))
end
fprintf('\n<strong>MS_Segment_raw</strong>: NLX csc and events loaded\n');

%% filtering Delta / theta

% delta filter.
% cfg_filt_d.display_filter = 1; % use this to see the fvtool and wait for a keypress before continuing.
delta_csc = FilterLFP(cfg.filt_d,csc);

% filter into the theta band
theta_csc = FilterLFP(cfg.filt_t, csc);


% get the theta_delta ratio
lfp_idx = find(strcmp(csc.label, 'LFP')); % use the LFP channel Or some other identifier in the csc.label. ATM this is to avoid using the 'EMG' channel.
emg_idx = find(strcmp(csc.label, 'EMG')); % use the LFP channel Or some other identifier in the csc.label. ATM this is to avoid using the 'EMG' channel.

td_ratio = abs(hilbert(theta_csc.data(lfp_idx, :))) ./ abs(hilbert(delta_csc.data(lfp_idx, :)));
td_ratio = filter(gausswin(4000)/sum(gausswin(4000)), 1, td_ratio); % 1d gaussian window with .5s window.

% narrow theta / wide theta (from Watson et al. 2016) [Note: that paper was
% in the neocortex.
t_emg_ratio = abs(hilbert(theta_csc.data(lfp_idx, :))) ./ abs(hilbert(csc.data(emg_idx, :)));
t_emg_ratio = filter(gausswin(4000)/sum(gausswin(4000)), 1, t_emg_ratio); % 1d gaussian window with .5s window.



% conv2(td_ratio,gausswin(1000),'same')
% conv2(td_ratio,gausskernel(1000,20),'same')

% add delta to the csc as a channel.
csc.data = cat(1,csc.data ,delta_csc.data(lfp_idx,:));
csc.label{end+1} = 'Delta';
%add theta
csc.data = cat(1,csc.data ,theta_csc.data(lfp_idx,:));
csc.label{end+1} = 'Theta';

% add theta-delta
csc.data = cat(1,csc.data ,td_ratio);
csc.label{end+1} = 'Theta/delta';

% add theta-wide theta
csc.data = cat(1,csc.data ,t_emg_ratio);
csc.label{end+1} = 'theta/emg';

clear delta_csc theta_csc td_ratio t_emg_ratio
fprintf('\n<strong>MS_Segment_raw</strong>: Delta, Theta, Theta./Delta ratio, and Theta./EMG have been added as csc channels.\n');



%% need to add a piece that will identify periods where the MS was recording but the NLX was not (example: when the mouse is on the track)
if length(evt_blocks) < length(TS)
    fprintf('Length of TS (%d) and evt_blocks (%d) are not equal. Checking for odd MS timestamp lengths out...\n',length(TS), length(evt_blocks));
    
    for iT = length(TS):-1:1
        l_ts(iT) = length(TS{iT}.system_clock{1});
    end
    
    for iE = length(evt_blocks):-1:1
        l_evt(iE) = length(evt_blocks{iE}.t{cfg.evt.t_chan});
    end
    
    odd_idx = find(~ismembertol(l_ts, l_evt, 5,'OutputAllIndices',true,'DataScale', 1));
    
    for iOdd = 1:length(odd_idx)
        fprintf('Found odd TS files at idx: %d.   Length: %d samples\n', odd_idx(iOdd), length(TS{odd_idx(iOdd)}.system_clock{1}));
    end
    
    % keep the track ms struct and save
    keep_idx = 1:size(ms_seg.RawTraces,1);
    keep_idx =keep_idx(find((keep_idx ~= odd_idx)));
    
    % get the pre, trk, post values
    for iC = length(TS):-1:1
        if iC < odd_idx
            ms_seg.pre_post{iC} = 'pre';
        elseif iC == odd_idx
            ms_seg.pre_post{iC} = 'task';
        elseif iC > odd_idx
            ms_seg.pre_post{iC} = 'post';
        end
    end
    
    cfg_rem = [];
    ms_trk = MS_remove_data_sandbox(cfg_rem, ms_seg, keep_idx);
    
    ms_trk = MS_de_cell(ms_trk);
    
    % binarize the trace
    
    ms_trk = msExtractBinary_detrendTraces(ms_trk);
    
    save([ms_resize_dir filesep 'ms_trk.mat'], 'ms_trk', '-v7.3')
    clear 'ms_trk';
    
    % remove from main ms struct
    cfg_rem = [];
    ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg, odd_idx);
    
    % remove from TS and labeling structs
    TS(odd_idx) = [];
    
    TS_name(odd_idx)= [];
    
    hypno_labels(odd_idx) = [];
    
    time_labels(odd_idx) = [];
    
    
    fprintf('\n<strong>MS_Segment_raw</strong>: miniscope epoch: %d was flagged for removal\n', odd_idx);
    
end
if ~isfield(ms_seg, 'removed')
    ms_seg.removed = {};
end
ms_seg.removed{end+1} = TS_name{odd_idx};


%% append the NLX data to the ms structure (be saure to use the same channel as the one used for extraction (cfg_evt_blocks.t_chan).
flag = [];
cut_offs = NaN(2,length(TS));
res_csc = cell(1, length(TS));
res_evt = cell(1,length(TS));

for iT = 1:length(TS)
    if length(TS{iT}.system_clock{1}) == length(evt_blocks{iT}.t{cfg.evt.t_chan})
        disp(['TS' num2str(iT) '-' TS_name{iT} ': ' num2str(length(TS{iT}.system_clock{1}))   'samples, '  num2str(length(TS{iT}.system_clock{1}) / TS{iT}.cfg.Fs{1},3) 'sec at ' num2str(TS{iT}.cfg.Fs{1},3) 'Hz'...
            'NLX: ' num2str(length(evt_blocks{iT}.t{cfg.evt.t_chan})) ' samples,' num2str(evt_duration(iT),3) 'at ' num2str(1/mode(diff(evt_blocks{iT}.t{cfg.evt.t_chan})),3) 'Hz'])
        res_csc{iT} = restrict(csc, evt_blocks{iT}.t{cfg.evt.t_chan}(1), evt_blocks{iT}.t{cfg.evt.t_chan}(end));
        res_evt{iT} = restrict(nlx_evts, evt_blocks{iT}.t{cfg.evt.t_chan}(1), evt_blocks{iT}.t{cfg.evt.t_chan}(end));
        
    elseif length(TS{iT}.system_clock{1}) - length(evt_blocks{iT}.t{cfg.evt.t_chan}) ==cfg.TS_nlx_match && cfg.TS_nlx_match ~= 0
        warning('TS do not match nlx .nev data. TS# %s  %s samples  - NLX: %s events',...
            num2str(iT), num2str(length(TS{iT}.system_clock{1})), num2str(length(evt_blocks{iT}.t{cfg.evt.t_chan})))
        fprintf('Times differ by %d sample(s). cfg.TS_nlx_match evoked, correcting TS...\n',abs(length(TS{iT}.system_clock{1}) - length(evt_blocks{iT}.t{cfg.evt.t_chan})))
        %         figure(121)
        %         subplot(2,1,1)
        %         plot(diff(TS{iT}.system_clock{1}))
        %         title(['TS timestamps. Segment: ' num2str(iT)])
        %         subplot(2,1,2)
        %         plot(diff(evt_blocks{iT}.t{cfg.evt.t_chan}))
        %         title('NLX evt timestamps')
        %         pause(1)
        
        %save values for resizing the ms struct data fields.
        cut_offs(:,iT) = [2;length(TS{iT}.system_clock{1})];
        % keep the nlx segment
        res_csc{iT} = restrict(csc, evt_blocks{iT}.t{cfg.evt.t_chan}(1), evt_blocks{iT}.t{cfg.evt.t_chan}(end));
        res_evt{iT} = restrict(nlx_evts, evt_blocks{iT}.t{cfg.evt.t_chan}(1), evt_blocks{iT}.t{cfg.evt.t_chan}(end));
        
    else
        warning('TS do not match nlx .nev data. TS# %s  %s samples  - NLX: %s events',...
            num2str(iT), num2str(length(TS{iT}.system_clock{1})), num2str(length(evt_blocks{iT}.t{cfg.evt.t_chan})))
        flag = [flag, iT];
        res_csc{iT} = [];
        res_evt{iT} = [];
    end
end

% remove unused blocks. In this case it is any one that does not match in
% the NLX events and the TS timestamps.
res_csc = res_csc(~cellfun('isempty',res_csc));
res_evt = res_evt(~cellfun('isempty',res_evt));


hypno_labels(flag) = [];
time_labels(flag) = [];
hypno_labels = hypno_labels(~cellfun('isempty', hypno_labels));
time_labels = time_labels(~cellfun('isempty', time_labels));
%% resize the ms_data to better match the NLX [soft fix for 1 sample differences]
if cfg.TS_nlx_match ==1
    cfg_resize = [];
    cfg_resize.tvec_to_use = 'time'; % could be 'time', or 'NLX_csc'
    cfg_resize.cutoffs = cut_offs; % should be [2 x nSegments] row 1 is start and row 2 is stop
    
    
    ms_seg = MS_resize_segments(cfg_resize, ms_seg);
    fprintf('\n<strong>MS_Segment_raw</strong>: cfg.TS_nlx_match evoked therefore: Resized ms segments')
    for iC = 1:length(cut_offs)
        if sum(isnan(cut_offs(1,iC)) + isnan(cut_offs(2,iC))) ~=2
            fprintf('<strong> %d</strong>', iC)
        end
        
    end
    
    fprintf(' that were 1 sample longer than the NLX seg.\n');
    
end


%% update the ms structure with the NLX data
cfg_rem = [];
% cfg_rem.user_fields = {'BinaryTraces'};
ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg, [flag]);
fprintf('\n<strong>MS_Segment_raw</strong>: miniscope epoch: %d was flagged for removal\n', flag);
for iR = 1:length(flag)
    ms_seg.removed{end+1} = TS_name{flag(iR)};
end
% add in the NLX data

ms_seg = MS_append_data_sandbox(ms_seg, 'NLX_csc', res_csc, 'NLX_evt', res_evt, 'hypno_label', hypno_labels, 'time_labels', time_labels);
fprintf('\n<strong>MS_Segment_raw</strong>: NLX_csc appended\n');

% clear large variables from workspace for memory.
% clear ms res_csc res_evt flag



%% get some emg stats for scaling
emg_chan  = find(ismember(cfg.csc.label, 'EMG')); % used to get the emg range.
% get the min and max emg range for the first 5mins of the recording. used for consistency.
cfg.resize = [];
cfg.resize.emg_range = [min(csc.data(emg_chan,1:(300*csc.cfg.hdr{1}.SamplingFrequency))), max(csc.data(emg_chan,1:(300*csc.cfg.hdr{1}.SamplingFrequency)))]; % get the min and max emg range for the first 10s of the recording. used for consistency.

%% spectrogram of an episode w/ability to resize using gui

[cut_vals, remove_flag] = MS_plot_spec_resize(cfg.resize, ms_seg);

flag = find(remove_flag);  % makes appending easier ina few steps
%% resize the events [WIP: has trouble resizing across ms and NLX timescales]
cfg_resize = [];
cfg_resize.tvec_to_use = 'NLX_csc'; % could be 'time', or 'NLX_csc'
cfg_resize.cutoffs = cut_vals; % should be [2 x nSegments] row 1 is start and row 2 is stop


ms_seg_resize = MS_resize_segments(cfg_resize, ms_seg);

%% remove segments that the user flagged in MS_plot_spec_resize
cfg_rem = [];
% cfg_rem.user_fields = {'BinaryTraces'};
ms_seg_resize = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, flag);
fprintf('\n<strong>MS_Segment_raw</strong>: miniscope epoch: %d was flagged for removal\n', find(remove_flag));

if ~isempty(flag)
    for iR = 1:length(flag)
        ms_seg.removed{end+1} = TS_name{flag(iR)};
    end
end



%% spectrogram of an episode w/
cfg.resize.resize = 0; % don't resize this time just plot.
MS_plot_spec_resize(cfg.resize, ms_seg_resize);


%% quick check?
%
% if isfield(cfg, 'check')
% %     cfg.check = [];
%     %     cfg_check.x_zoom = [ 0 5];
% %     cfg.check.Ca_type = 'RawTraces';
%     %     cfg_check.Ca_type = 'BinaryTraces';
%     cfg.check.chan_to_plot = ms_seg.NLX_csc{1}.label;
%     cfg.check.plot_type = '3d';
%     cfg.check.label = 'hypno_label';
%     cfg.check.emg_range = emg_range;
%
%     MS_plot_ca_nlx(cfg.check, ms_seg, res_csc);
% end

%% binarize the traces in each segment and save each one back to the same folder name as the original Ms TS file.
% set up empty variables for each PRE v POST and REM v SW
all_binary_pre = []; all_binary_post= [];
all_RawTraces_pre = []; all_RawTraces_post = [];
all_detrendRaw_pre = []; all_detrendRaw_post = [];
all_seg_idx = [];
%rem
all_binary_pre_REM = []; all_binary_post_REM= [];
all_RawTraces_pre_REM = []; all_RawTraces_post_REM = [];
all_detrendRaw_pre_REM = []; all_detrendRaw_post_REM = [];
pre_REM_idx = []; post_REM_idx = [];
% SW
all_binary_pre_SW = []; all_binary_post_SW= [];
all_RawTraces_pre_SW = []; all_RawTraces_post_SW = [];
all_detrendRaw_pre_SW = []; all_detrendRaw_post_SW = [];
pre_SW_idx = []; post_SW_idx = [];


for iSeg = 1:length(ms_seg_resize.RawTraces)
    
    keep_idx = 1:size(ms_seg_resize.RawTraces,1);
    keep_idx =keep_idx(find((keep_idx ~= iSeg)));
    
    cfg_rem = [];
    this_ms = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
    
    this_ms = MS_de_cell(this_ms);
    
    % binarize the trace
    
    this_ms = msExtractBinary_detrendTraces(this_ms);
    
    this_dir = [];
    this_dir = [ms_resize_dir filesep ms_seg_resize.file_names{iSeg}];
    fprintf('<strong>%s</strong>: saving resized ms struct back to %s...\n', mfilename, this_dir)
    mkdir(this_dir)
    save([this_dir filesep 'ms_seg_resize_' ms_seg_resize.pre_post{iSeg} '_' ms_seg_resize.hypno_label{iSeg}], 'this_ms', '-v7.3');
    
    % keep the index for the segment.
    if isempty(all_seg_idx)
        all_seg_idx(iSeg) = length(this_ms.RawTraces);
    else
        all_seg_idx(iSeg) = length(this_ms.RawTraces) + all_seg_idx(iSeg -1);
    end
    % cat the binary traces for pre V post, and REM v SW
    if strcmp(ms_seg_resize.pre_post{iSeg}, 'pre')
        all_binary_pre = [all_binary_pre; this_ms.Binary];
        all_RawTraces_pre = [all_RawTraces_pre; this_ms.RawTraces];
        all_detrendRaw_pre = [all_detrendRaw_pre; this_ms.detrendRaw];
        
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            all_binary_pre_REM = [all_binary_pre_REM; this_ms.Binary];
            all_RawTraces_pre_REM = [all_RawTraces_pre_REM; this_ms.RawTraces];
            all_detrendRaw_pre_REM = [all_detrendRaw_pre_REM; this_ms.detrendRaw];
            pre_REM_idx = [pre_REM_idx, iSeg];
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            all_binary_pre_SW = [all_binary_pre_SW; this_ms.Binary];
            all_RawTraces_pre_SW = [all_RawTraces_pre_SW; this_ms.RawTraces];
            all_detrendRaw_pre_SW = [all_detrendRaw_pre_SW; this_ms.detrendRaw];
            pre_SW_idx = [pre_SW_idx, iSeg];
        end
        
    elseif strcmp(ms_seg_resize.pre_post{iSeg}, 'post')
        all_binary_post = [all_binary_post; this_ms.Binary];
        all_RawTraces_post = [all_RawTraces_post; this_ms.RawTraces];
        all_detrendRaw_post = [all_detrendRaw_post; this_ms.detrendRaw];
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            all_binary_post_REM = [all_binary_post_REM; this_ms.Binary];
            all_RawTraces_post_REM = [all_RawTraces_post_REM; this_ms.RawTraces];
            all_detrendRaw_post_REM = [all_detrendRaw_post_REM; this_ms.detrendRaw];
            post_REM_idx = [post_REM_idx, iSeg];
            
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            all_binary_post_SW = [all_binary_post_SW; this_ms.Binary];
            all_RawTraces_post_SW = [all_RawTraces_post_SW; this_ms.RawTraces];
            all_detrendRaw_post_SW = [all_detrendRaw_post_SW; this_ms.detrendRaw];
            post_SW_idx = [post_SW_idx, iSeg];
            
        end
    end
end
all_seg_idx = [0 all_seg_idx];

fprintf('<strong>%s</strong>: saving concatinating Binary, RawTraces, detrendRaw, and indicies', mfilename);

% save everything
save([ms_resize_dir filesep 'all_seg_idx.mat'], 'all_seg_idx', '-v7.3');

% pre only
save([ms_resize_dir filesep 'all_binary_pre.mat'], 'all_binary_pre', '-v7.3');
save([ms_resize_dir filesep 'all_RawTraces_pre.mat'], 'all_RawTraces_pre', '-v7.3');
save([ms_resize_dir filesep 'all_detrendRaw_pre.mat'], 'all_detrendRaw_pre', '-v7.3');

% post only
save([ms_resize_dir filesep 'all_binary_post.mat' ], 'all_binary_post', '-v7.3');
save([ms_resize_dir filesep 'all_RawTraces_post.mat'], 'all_RawTraces_post', '-v7.3');
save([ms_resize_dir filesep 'all_detrendRaw_post.mat'], 'all_detrendRaw_post', '-v7.3');


% pre REM only
save([ms_resize_dir filesep 'all_binary_pre_REM.mat'], 'all_binary_pre_REM', '-v7.3');
save([ms_resize_dir filesep 'all_RawTraces_pre_REM.mat'], 'all_RawTraces_pre_REM', '-v7.3');
save([ms_resize_dir filesep 'all_detrendRaw_pre_REM.mat'], 'all_detrendRaw_pre_REM', '-v7.3');

% post SW only
save([ms_resize_dir filesep 'all_binary_post_REM.mat' ], 'all_binary_post_REM', '-v7.3');
save([ms_resize_dir filesep 'all_RawTraces_post_REM.mat'], 'all_RawTraces_post_REM', '-v7.3');
save([ms_resize_dir filesep 'all_detrendRaw_post_REM.mat'], 'all_detrendRaw_post_REM', '-v7.3');

% pre REM only
save([ms_resize_dir filesep 'all_binary_pre_SW.mat'], 'all_binary_pre_SW', '-v7.3');
save([ms_resize_dir filesep 'all_RawTraces_pre_SW.mat'], 'all_RawTraces_pre_SW', '-v7.3');
save([ms_resize_dir filesep 'all_detrendRaw_pre_SW.mat'], 'all_detrendRaw_pre_SW', '-v7.3');

% post SW only
save([ms_resize_dir filesep 'all_binary_post_SW.mat' ], 'all_binary_post_SW', '-v7.3');
save([ms_resize_dir filesep 'all_RawTraces_post_SW.mat'], 'all_RawTraces_post_SW', '-v7.3');
save([ms_resize_dir filesep 'all_detrendRaw_post_SW.mat'], 'all_detrendRaw_post_SW', '-v7.3');

% visualize
figure(1010)
subplot(3,1,1)
title('Binary Pre cat: SW: green REM: red')
hold on
for ii = 1:10
    plot(all_binary_pre(:,ii)+ii);
end
vline(all_seg_idx(pre_REM_idx), {'r'});
vline(all_seg_idx(pre_SW_idx),{'g'})

subplot(3,1,2)
title('RawTraces Pre cat: SW: green REM: red')
hold on
for ii = 1:10
    plot(all_RawTraces_pre(:,ii)+ii);
end
vline(all_seg_idx(pre_REM_idx), {'r'});
vline(all_seg_idx(pre_SW_idx),{'g'})

subplot(3,1,3)
title('detrendRaw Pre cat: SW: green REM: red')

hold on
for ii = 1:10
    plot(all_detrendRaw_pre(:,ii)+ii);
end
vline(all_seg_idx(pre_REM_idx), {'r'});
vline(all_seg_idx(pre_SW_idx),{'g'})
pause(5)
close; 
%% clean up and export the ms_seg_resize
save([ms_resize_dir filesep 'ms_resize.mat'], 'ms_seg_resize', '-v7.3')