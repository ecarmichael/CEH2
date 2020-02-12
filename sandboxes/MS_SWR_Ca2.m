% SWR detection and clustering sandbox

%% add paths

close all
restoredefaultpath
global PARAMS
os = computer;

if ismac
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.ft_code_dir = '/Users/jericcarmichael/Documents/GitHub/fieldtrip'; % FieldTrip toolbbox (used for spectrogram visualization)
    
elseif strcmp(os, 'GLNXA64')
    
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/'; % raw data location.
    
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.chronux_code_dir = '/home/ecarmichael/Documents/chronux/chronux_2_12'; % FieldTrip toolbbox (used for spectrogram visualization)
    PARAMS.ft_code_dir = '/home/ecarmichael/Documents/GitHub/fieldtrip'; % FieldTrip toolbbox (used for spectrogram visualization)
    PARAMS.seqNMF_dir = '/home/ecarmichael/Documents/GitHub/seqNMF';% seqNMF pathway. used for sequence detection.
    
else
    disp('on a PC fill this in yourself....')
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

% try the newer NLX loaders for UNIX
[~, d] = version;
if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os


%% load the Miniscope data

load('ms.mat')

% get the timestampsaddpath(PARAMS.seqNMF_dir)

raw_data_folder = strsplit(PARAMS.data_dir, filesep);
[TS, TS_name] = MS_collect_timestamps(strjoin([PARAMS.raw_data_dir raw_data_folder(end)], ''));

% [TS, TS_name] = MS_collect_timestamps('/home/ecarmichael/Documents/Williams_Lab/7_12_2019_PV1069_LTD5');

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
%% Get some basic statistics for the individual Ca Traces
cfg_stats = [];
cfg_stats.data_type = 'RawTraces';

ms = MS_characterize_trace(cfg_stats, ms);
fprintf('\n<strong>MS_SWR_Ca2</strong>: basic statstics computed for each Ca trace');


%% remove cells with low quality signals.
cfg_remove_trace = [];
cfg_remove_trace.threshold = 8;
cfg_remove_trace.remove_idx = ms.stats.RawTraces.z_max <= cfg_remove_trace.threshold;

ms = MS_Remove_trace(cfg_remove_trace, ms);
fprintf('\n<strong>MS_SWR_Ca2</strong>: %d Traces had max zscore values  < %d. Taken to indicate poor spike quality \n',sum(cfg_remove_trace.remove_idx), cfg_remove_trace.threshold);



%%  convert the Ca transitents into a binarized vector
cfg_bin = [];
cfg_bin.method = 'rise';
cfg_bin.rise.smooth_type = 'sgolay';
cfg_bin.threshold = 2;
ms = MS_binarize_data_sandbox(cfg_bin, ms);
fprintf('\n<strong>MS_SWR_Ca2</strong>: miniscope data has been binarized using a %s method with a threshold of %d\n', cfg_bin.method, cfg_bin.threshold);
pause(2)
close all


%% segment the data
cfg_seg = [];
cfg_seg.user_fields = {'BinaryTraces'};
ms_seg = MS_segment_ms_sandbox(cfg_seg, ms);

fprintf('\n<strong>MS_SWR_Ca2</strong>: miniscope data has been segmented into %d individual recording epochs\n method used: %s\n', length(ms_seg.time), ms_seg.format);

%% Load nlx data

% load the Keys file with all of the experiment details.
%(can be generated with the 'MS_Write_Keys' function)
if exist('*Keys.m', 'file') 
    ExpKeys = MS_Load_Keys();
end

% load events
nlx_evts = LoadEvents([]); % get '.nev' file in this dir.

% load the NLX CSC data (using vandermeer lab code) [todo:replace with own]
cfg_csc = [];
cfg_csc.fc = {'CSC7.ncs', 'CSC1.ncs'}; % use csc files from Keys. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
cfg_csc.desired_sampling_frequency = 2000; 
csc = MS_LoadCSC(cfg_csc); % need to comment out ExpKeys lines in LoadCSC
csc.label{1} = 'LFP'; 
csc.label{2} = 'EMG'; 

% extract NLX event epochs
cfg_evt_blocks = [];
cfg_evt_blocks.t_chan = 5;
nlx_evts.t{5} = sort([nlx_evts.t{3}, nlx_evts.t{4}]);
nlx_evts.label{5} = 'merge TTls at 3 and 4';
[evt_blocks, evt_iv, evt_duration] = MS_extract_NLX_blocks_sandbox(cfg_evt_blocks, nlx_evts);


% compare to TS to ms
fprintf('\n****Comparing TS files to processed miniscope (ms) data\n')
if length(TS) ~= length(evt_blocks)
    warning('Number of Timestamp files (%s) does not match the number of detected NLX event blocks (%s)', num2str(length(TS)),num2str(length(evt_blocks)))
end

%% need to add a piece that will identify periods where the MS was recording but the NLX was not (example: when the mouse is on the track)
if length(evt_blocks) < length(TS)
    fprintf('Length of TS (%d) and evt_blocks (%d) are not equal. Checking for odd MS timestamp lengths out...\n',length(TS), length(evt_blocks));
    
    for iT = length(TS):-1:1
        l_ts(iT) = length(TS{iT}.system_clock{1});
    end
    
    for iE = length(evt_blocks):-1:1
        l_evt(iE) = length(evt_blocks{iE}.t{cfg_evt_blocks.t_chan});
    end

    odd_idx = find(~ismembertol(l_ts, l_evt, 5,'OutputAllIndices',true,'DataScale', 1));
    
    for iOdd = 1:length(odd_idx)
     fprintf('Found odd TS files at idx: %d.   Length: %d samples\n', odd_idx(iOdd), length(TS{odd_idx(iOdd)}.system_clock{1}));
    end
    
    % remove from ms struct
    cfg_rem = [];
    cfg_rem.user_fields = {'BinaryTraces'};
    ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg, [odd_idx]);
    
    % remove from TS and labeling structs
    TS(odd_idx) = [];
%     TS = TS(~cellfun('isempty',TS));
    
    TS_name(odd_idx)= [];
%     TS_name = TS_name(~cellfun('isempty',TS_name));
    
    hypno_labels(odd_idx) = [];
%     hypno_labels = hypno_labels(~cellfun('isempty',hypno_labels));

    time_labels(odd_idx) = [];
%     time_labels = time_labels(~cellfun('isempty',time_labels));


    fprintf('\n<strong>MS_SWR_Ca2</strong>: miniscope epoch: %d was flagged for removal\n', odd_idx);

end


%% append the NLX data to the ms structure (be saure to use the same channel as the one used for extraction (cfg_evt_blocks.t_chan).
flag = [];
res_csc = cell(1, length(TS));
res_evt = cell(1,length(TS));
for iT = 1:length(TS)
    if length(TS{iT}.system_clock{1}) == length(evt_blocks{iT}.t{cfg_evt_blocks.t_chan})
        disp(['TS' num2str(iT) '-' TS_name{iT} ': ' num2str(length(TS{iT}.system_clock{1}))   'samples, '  num2str(length(TS{iT}.system_clock{1}) / TS{iT}.cfg.Fs{1},3) 'sec at ' num2str(TS{iT}.cfg.Fs{1},3) 'Hz'...
            'NLX: ' num2str(length(evt_blocks{iT}.t{cfg_evt_blocks.t_chan})) ' samples,' num2str(evt_duration(iT),3) 'at ' num2str(1/mode(diff(evt_blocks{iT}.t{cfg_evt_blocks.t_chan})),3) 'Hz'])
        res_csc{iT} = restrict(csc, evt_blocks{iT}.t{cfg_evt_blocks.t_chan}(1), evt_blocks{iT}.t{cfg_evt_blocks.t_chan}(end));
        res_evt{iT} = restrict(nlx_evts, evt_blocks{iT}.t{cfg_evt_blocks.t_chan}(1), evt_blocks{iT}.t{cfg_evt_blocks.t_chan}(end));
        
    else
        warning('TS do not match nlx .nev data. TS# %s  %s samples  - NLX: %s events',...
            num2str(iT), num2str(length(TS{iT}.system_clock{1})), num2str(length(evt_blocks{iT}.t{cfg_evt_blocks.t_chan})))
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


%% update the ms structure with the NLX data
cfg_rem = [];
cfg_rem.user_fields = {'BinaryTraces'};
ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg, [flag]);
fprintf('\n<strong>MS_SWR_Ca2</strong>: miniscope epoch: %d was flagged for removal\n', flag);

ms_seg = MS_append_data_sandbox(ms_seg, 'NLX_csc', res_csc, 'NLX_evt', res_evt, 'hypno_label', hypno_labels, 'time_labels', time_labels);
fprintf('\n<strong>MS_SWR_Ca2</strong>: NLX_csc appended\n');

% clear large variables from workspace for memory.
% clear ms res_csc res_evt flag




%% quick check?
check = 1; % toggle to skip check plots.

if check  ==1
    cfg_check = [];
    %     cfg_check.x_zoom = [ 0 5];
        cfg_check.Ca_type = 'RawTraces';
%     cfg_check.Ca_type = 'BinaryTraces';
    cfg_check.plot_type = '3d';
    cfg_check.label = 'hypno_label';
    MS_plot_ca_nlx(cfg_check, ms_seg, res_csc);
    
    
%     figure
%     subplot(6,6,1:5)
%     plot(ms.time/1000,smoothdata(sum(ms.BinaryTraces,2),'gaussian',1000), 'linewidth', 4)
%     hline(0, 'k')
%     color = get(gcf,'Color');
%     set(gca, 'color', color);
%     %     set(gca,'XColor',color)%,'YColor',color,'TickDir','out')
%     xlim([ms.time(1)/1000 ms.time(end)/1000]);
%     xticks([]);
%     title('Sum by time')
%     subplot(6,6,[12 18 24 30 36])
%     plot(smoothdata(sum(ms.BinaryTraces,1),'gaussian',25),1:size(ms.BinaryTraces,2), 'linewidth', 4)
%     yticks([]); ylim([1 size(ms.BinaryTraces,2)]);
%     ylabel('Sum by cells')
%     set(gca, 'color', color);
%     yyaxis right
%     hax = gca;
%     set(hax.YAxis, {'color'}, {'k'})
%     %     set(gca,'XColor',color,'YColor',color,'TickDir','out')
%     ylim([1 size(ms.BinaryTraces,2)]);
%     subplot(6, 6, [7 13 19 25 31 8 14 20 26 32 9 15 21 27 33 10 16 22 28 34 11 17 23 29 35])
%     imagesc(ms.time/1000,1:size(ms.BinaryTraces,2), ms.BinaryTraces)
%     colormap(flipud(hot))
%     ylabel('Cell #'); xlabel('Time (s)');
%     xlim([ms.time(1)/1000 ms.time(end)/1000]);
    
    
end

%% spectrogram of an episode

this_block = 10
[~,F,T,P] = spectrogram(ms_seg.NLX_csc{this_block}.data, rectwin(128), 128/4, 1:.1:40, csc.cfg.hdr{1}.SamplingFrequency);

figure(1112)
ax1 = imagesc(T,F,10*log10(P));
set(ax1, 'AlphaData', ~isinf(10*log10(P)))
set(gca,'FontSize',28);
axis xy xlabel('Time (s)'); ylabel('Frequency (Hz)');
ax = gca;
% ax.YTickLabels = {0 [], [], 60, [], [], 80}; % set
set(gca, 'tickdir', 'out')


hold on

plot(ms_seg.NLX_csc{this_block}.tvec-ms_seg.NLX_csc{this_block}.tvec(1) , (ms_seg.NLX_csc{this_block}.data*10000)+20, 'k' )


cfg_filt = [];
% cfg_filt.f = [5 11]; %setting theta (hertz)
cfg_filt.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt.f  = [1 5];
cfg_filt.order = 8; %type filter order
cfg_filt.display_filter = 1; % use this to see the fvtool
delta_csc = FilterLFP(cfg_filt,ms_seg.NLX_csc{this_block});

% filter into the theta band
cfg_filt = [];
% cfg_filt.f = [5 11]; %setting theta (hertz)
cfg_filt.type = 'cheby1'; %the type of filter I want to use via filterlfp
cfg_filt.f  = [6 11];
cfg_filt.order = 8; %type filter order
cfg_filt.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  Looks like trash with ord ~= 3 in cheby1. Butter and fdesign are also trash.
theta_csc = FilterLFP(cfg_filt, ms_seg.NLX_csc{this_block});


%% segment data into one of the specified recording blocks should be hard

SW_block = 4; % good block based on visual inspection of the check plots with hypno labels (above)
REM_block = 9; % nice REM block

fprintf('\n<strong>MS_SWR_Ca2</strong>: using recording blocks <strong>SW = %d (%.1fs) REM = %d (%.1fs)</strong>\n', SW_block,(ms_seg.time{SW_block}(end) - ms_seg.time{SW_block}(1))*0.001, REM_block, (ms_seg.time{REM_block}(end) - ms_seg.time{REM_block}(1))*0.001)

%% extract SWR candidate events

for iBlock = [SW_block]%, REM_block]
    
    this_csc = res_csc{iBlock}; % pull out a block of
    
    
    
    
    
    
    
    %% basic filtering and thresholding
    % mouse SWR parameters are based off of Liu, McAfee, & Heck 2017 https://www.nature.com/articles/s41598-017-09511-8#Sec6
    check = 1; % used for visual checks on detected events.
    ft_check = 1; % use fieldtrip to
    %set up ripple band
    cfg_filt = [];
    cfg_filt.type = 'butter'; %Cheby1 is sharper than butter
    cfg_filt.f  = [140 250]; % broad, could use 150-200?
    cfg_filt.order = 4; %type filter order (fine for this f range)
    cfg_filt.display_filter = 0; % use this to see the fvtool
    csc_ripple = FilterLFP(cfg_filt, this_csc);
    
    
    % convert to amplitude or power
    amp_ripple = csc_ripple; % clone to make things simple and replace
    
    
    for iChan = 1:size(csc_ripple.data,1)
        amp_ripple.data(iChan,:) = abs(hilbert(csc_ripple.data(iChan,:)));
        % Convolve with a gaussian kernel (improves detection)
        kernel = gausskernel(60,20); % note, units are in samples; for paper Methods, need to specify Gaussian SD in ms
        fprintf('\nGausskernal using 60 samples = %0.0fms with SD = 20 samples (%0.0fms)\n', (60/this_csc.cfg.hdr{1}.SamplingFrequency)*1000, (20/this_csc.cfg.hdr{1}.SamplingFrequency)*1000)
        amp_ripple.data(iChan,:) = conv(amp_ripple.data(iChan,:),kernel,'same');
        amp_ripple.units = 'amplitude';
        % if you want units in power use:
        %     amp_ripple.data(iChan,:) =  amp_ripple.data(iChan,:).^2;
        %     amp_ripple.units = 'power';
        
    end
    
    if check
        figure(111)
        plot(this_csc.tvec, this_csc.data(1,:),'k',csc_ripple.tvec, csc_ripple.data(1,:), 'r',...
            amp_ripple.tvec, amp_ripple.data(1,:),'b')
        legend({'Raw', '140-250 filt', 'Amp'})
        pause(1); close;
    end
    
    %% remove large amplitude artifacts before SWR detection
    
    
    csc_artif = this_csc;
    for iChan = 1:size(csc_ripple.data,1)
        csc_artif.data(iChan,:) = abs(csc_artif.data(iChan,:)); % detect artifacts both ways
    end
    
    cfg_artif_det = [];
    cfg_artif_det.method = 'raw';
    cfg_artif_det.threshold = std(csc_artif.data(1,:))*5;
    % cfg_artif_det.minlen = 0.01;
    cfg_artif_det.target = this_csc.label{1};
    evt_artif = TSDtoIV(cfg_artif_det,csc_artif);
    
    cfg_temp = []; cfg_temp.d = [-0.5 0.5];
    artif_evts = ResizeIV(cfg_temp,evt_artif);
    
    
    % plot
    if check
        plot(113)
        cfg_plot=[];
        cfg_plot.display = 'iv'; % tsd, iv
        cfg_plot.target = this_csc.label{1};
        PlotTSDfromIV(cfg_plot,artif_evts,csc_artif);
        hline(cfg_artif_det.threshold )
        pause(3); close all;
    end
    
    % zero pad artifacts to improve reliability of subsequent z-scoring
    artif_idx = TSD_getidx2(this_csc,evt_artif); % if error, try TSD_getidx (slower)
    for iChan = 1:size(csc_ripple.data,1)
        csc_ripple.data(iChan,artif_idx) = 0;
        amp_ripple.data(iChan,artif_idx) = 0;
    end
    
    % plot
    if check
        %     plot(114)
        hold on
        plot(amp_ripple.tvec, csc_ripple.data(1,:),'g');
        plot(amp_ripple.tvec, amp_ripple.data(1,:),'-k');
        
        pause(3); close all;
    end
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d large amplitude artifacts detected and zero-padded from csc_ripple.\n',length(artif_evts.tstart));
    
    
    %% isolate candidate events
    
    % get the thresholds
    cfg_detect = [];
    cfg_detect.operation = '>';
    cfg_detect.dcn = cfg_detect.operation; % b/c odd var naming in TSDtoIV
    cfg_detect.method = 'zscore';
    cfg_detect.threshold = 2;
    cfg_detect.target = csc.label{1};
    cfg_detect.minlen = 0.020; % 40ms from Vandecasteele et al. 2015
    cfg_detect.merge_thr = 0.02; % merge events that are within 20ms of each other.
    
    [swr_evts,evt_thr] = TSDtoIV(cfg_detect,amp_ripple);
    
    % % now apply to all data
    % cfg_select = [];
    % cfg_select.dcn = '>';
    % cfg_select.method = 'raw';
    % cfg_select.threshold = evt_thr;
    % cfg_select.target = 'CSC1.ncs';
    % cfg_select.minlen = cfg_detect.minlen;
    %
    % [evt_ids,~] = TSDtoIV(cfg_select,amp_ripple);
    
    
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events detected initially.\n',length(swr_evts.tstart));
    
    if check
        cfg_plot = [];
        cfg_plot.display = 'iv';
        cfg_plot.mode = 'center';
        cfg_plot.width = 0.2;
        cfg_plot.target = csc.label{1};
        
        PlotTSDfromIV(cfg_plot,swr_evts,csc);
        pause(2); close all;
    end
    
    
    
    %% exclude events with insufficient cycles - count how many exist above same threshold as used for detection
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filt;
    swr_evt_out = CountCycles(cfg_cc,this_csc,swr_evts);
    
    % get get the evetns with sufficient cycles.
    cfg_gc = [];
    cfg_gc.operation = '>=';
    cfg_gc.threshold = 4;
    swr_evt_out = SelectIV(cfg_gc,swr_evt_out,'nCycles');
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events remain after cycle count thresholding (%d cycle minimum).\n',length(swr_evt_out.tstart), cfg_gc.threshold);
    
    %% check for evnts that are too long.
    % add in a user field for the length of the events (currently not used)
    swr_evt_out.usr.evt_len = (swr_evt_out.tend - swr_evt_out.tstart)';
    
    cfg_max_len = [];
    cfg_max_len.operation = '<';
    cfg_max_len.threshold = .1;
    swr_evt_out = SelectIV(cfg_max_len,swr_evt_out,'evt_len');
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>:: %d events remain after event length cutoff (> %d ms removed).\n',length(swr_evt_out.tstart), (cfg_max_len.threshold)*1000);
    
    
    %% check for evnts with high raw varience. 'var_raw' is added as a swr_evt_out.usr field in CountCycles
    
    cfg_max_len = [];
    cfg_max_len.operation = '<';
    cfg_max_len.threshold = 1;
    swr_evt_out = SelectIV(cfg_max_len,swr_evt_out,'var_raw');
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events remain after raw varience thresholding (''var_raw'' > %d removed).\n',length(swr_evt_out.tstart), cfg_max_len.threshold);
    
    %% remove events that cooinside with artifacts.
    swr_evt_out = DifferenceIV([], swr_evt_out, artif_evts);
    
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %d events remain after removing those co-occuring with artifacts.\n',length(swr_evt_out.tstart));
    
    %% check again
    if check
        cfg_plot = [];
        cfg_plot.display = 'iv';
        cfg_plot.mode = 'center';
        cfg_plot.width = 0.2;
        cfg_plot.target = csc.label{1};
        cfg_plot.title = 'var';
        PlotTSDfromIV(cfg_plot,swr_evt_out,csc);
        pause(3); close all;
    end
    
    
    
    
    %% using FieldTrip Toolbox  (https://github.com/fieldtrip)
    %
    if ft_check == 1
        addpath(PARAMS.ft_code_dir);
        
        ft_defaults
        %
        % fc = {'CSC7.ncs'};
        % data = ft_read_neuralynx_interp(fc); used to updae TSDtoFT to give
        % correct formating. Works as MS_TSDtoFT.
        
        % convert data to ft format and turn into trials.
        data_ft = MS_TSDtoFT([], this_csc); % convert to ft format.
        
        swr_centers = IVcenters(swr_evt_out); % get the center of the swr event.
        
        cfg_trl = [];
        cfg_trl.t = cat(1,swr_centers);
        cfg_trl.t = cfg_trl.t - data_ft.hdr.FirstTimeStamp;
        cfg_trl.twin = [-1 1];
        cfg_trl.hdr = data_ft.hdr;
        
        trl = ft_maketrl(cfg_trl);
        
        cfg = [];
        cfg.trl = trl;
        data_trl = ft_redefinetrial(cfg,data_ft);
        
        
        cfg              = []; % start with empty cfg
        cfg.output       = 'pow';
        cfg.channel      = data_ft.label{1};
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 20:2:250; % frequencies of interest
        cfg.t_ftimwin    = ones(size(cfg.foi)).*0.05;%20./cfg.foi;  % window size: fixed at 0.5s
        cfg.toi          = -.2:0.0025:0.2; % times of interest
        cfg.pad          = 'nextpow2'; % recommened by FT to make FFT more efficient.
        
        TFR = ft_freqanalysis(cfg, data_trl);
        
        % track config for plotting.
        freq_params_str = sprintf('Spec using %0.0d swrs. Method: %s, Taper: %s', length(trl),cfg.method, cfg.taper);
        
        figure(iBlock)
        subplot(2,3,[1 2 4 5])
        cfg = [];
        cfg.channel      = data_ft.label{1};
        cfg.baseline     = [-1 -.01];
        cfg.baselinetype = 'relative';
        cfg.title = freq_params_str;
        ft_singleplotTFR(cfg, TFR);
        
        
        clear data_ft data_trl swr_centers
    end % end of ft_check
    
    

    
    for iBand = [3 6 ]%9]
        if iBand ==6
            data_ft = MS_TSDtoFT([], this_csc); % convert to ft format.
            title_name = 'Raw';
%             these_times = swr_evt_out.tstart;

        elseif iBand ==3

            data_ft = MS_TSDtoFT([], csc_ripple); % convert to ft format.
            title_name = 'Ripple band';
            % get the peak of the SWR for alignment
            for iTrial = length(swr_evt_out.tstart):-1:1
                this_swr = restrict(csc_ripple, swr_evt_out);
                [~, pk_idx] = max(this_swr.data);
                SWR_peak_t(iTrial) = this_swr.tvec(pk_idx); 
                these_times = SWR_peak_t;
            end
        end
        
            % average SWR LFP
            cfg_trl = [];
            %             cfg_trl.t = swr_evt_out.tstart;
            cfg_trl.t = these_times;
            cfg_trl.t = cfg_trl.t - data_ft.hdr.FirstTimeStamp;
            cfg_trl.twin = [-0.25 0.25];
            cfg_trl.hdr = data_ft.hdr;
            
            trl = ft_maketrl(cfg_trl);
            cfg = [];
            cfg.trl = trl;
            data_trl = ft_redefinetrial(cfg,data_ft);
        
  
        for iTrial = length(data_trl.trial):-1:1
            %     plot(data_trl.time{iTrial}, data_trl.trial{iTrial})
            disp(num2str(length(data_trl.trial{iTrial})));
            if length(data_trl.trial{iTrial}) < length(all_trials) % nan_pad if different which only seemed to happen once.
                pad_len = length(all_trials) - length(data_trl.trial{iTrial});
                all_trials(iTrial,:) = [data_trl.trial{iTrial}, NaN(pad_len,1)];
            else
                
                all_trials(iTrial,:) = data_trl.trial{iTrial};
            end
            %     pause(.5)
        end
        subplot(2,3,iBand)
        shadedErrorBar(data_trl.time{1}, nanmean(all_trials,1), nanstd(all_trials,1))
        title(title_name)
    end
    
            rmpath(PARAMS.ft_code_dir);

    %% block clean up
    
    
    % summary
    fprintf('\n<strong>MS_SWR_Ca2</strong>: %.0f candidate events.  Rate: %.1f/min , mean duration: %.1fms\n',length(swr_evt_out.tstart), length(swr_evt_out.tstart)/(((ms_seg.time{SW_block}(end) - ms_seg.time{SW_block}(1))*0.001)/60), (mean([swr_evt_out.tend-swr_evt_out.tstart]))*1000);
    
    
    SWR_candidates.(ms_seg.hypno_label{iBlock}) = swr_evt_out;
    
    
    clear this_csc artif_evts artif_idx csc_ripple csc_artif
    
end % end recording block iBlock

%% Ca coactivity with SWRs

swr_centers = IVcenters(SWR_candidates.SW); % convert to centered events;

% get the time idx that matches the SWR centers (use this if you just want
% one frame before and one after or something.
swr_ms_idx_centers = nearest_idx3(swr_centers, ms_seg.NLX_evt{SW_block}.t{cfg_evt_blocks.t_chan});

% alternative:
% get the idx for the start and end of the event.
swr_ms_idx_tstart = nearest_idx3(SWR_candidates.SW.tstart, ms_seg.NLX_evt{SW_block}.t{cfg_evt_blocks.t_chan});

swr_ms_idx_tend = nearest_idx3(SWR_candidates.SW.tend, ms_seg.NLX_evt{SW_block}.t{cfg_evt_blocks.t_chan});

% initialize some matricies to store the co-activity.
co_mat = NaN(size(ms_seg.BinaryTraces{SW_block},2),size(ms_seg.BinaryTraces{SW_block},2),length(swr_ms_idx_tstart)); % Make an empty matrix for co-activity
corr_mat = NaN(size(ms_seg.BinaryTraces{SW_block},2),size(ms_seg.BinaryTraces{SW_block},2),length(swr_ms_idx_tstart)); % Make an empty matrix for correlations
SWR_activity = NaN(size(ms_seg.BinaryTraces{SW_block},2),length(swr_ms_idx_tstart));

SWR_cat_data = [];
idx_win = [-50 50]; % window (in index values) around the event.

% initialize an image to draw over
figure(111)
h = imagesc(co_mat(:,:,1));

for iE =1: length(swr_ms_idx_tstart)%:-1:1 % loop SWRS
    
    if swr_ms_idx_tstart(iE)+idx_win(1) <= 0
        this_evt = ms_seg.BinaryTraces{SW_block}(1:swr_ms_idx_tstart(iE)+idx_win(2),:); % get all the values within this event.
        
    elseif swr_ms_idx_tstart(iE)+idx_win(2) >= length(ms_seg.BinaryTraces{SW_block})
        this_evt = ms_seg.BinaryTraces{SW_block}(swr_ms_idx_tstart(iE)+idx_win(1):end,:); % get all the values within this event.
        
    else
        this_evt = ms_seg.BinaryTraces{SW_block}(swr_ms_idx_tstart(iE)+idx_win(1):swr_ms_idx_tstart(iE)+idx_win(2),:); % get all the values within this event.
    end
    SWR_activity(:,iE) = sum(this_evt,1) >0; % see if anything was active.
    
    %     corr_mat(:,:,iE) = corr(SWR_activity', 'rows', 'pairwise');
    
    for ii = length(SWR_activity(:,iE)):-1:1
        for jj = length(SWR_activity(:,iE)):-1:1
            if SWR_activity(ii,iE) == 1 && SWR_activity(jj,iE) == 1
                co_mat(ii,jj,iE) = 1;
            elseif (SWR_activity(ii,iE) == 1 && SWR_activity(jj,iE) == 0) || (SWR_activity(ii,iE) == 0 && SWR_activity(jj,iE) == 1)
                co_mat(ii,jj,iE) = 0;
            elseif SWR_activity(ii,iE) == 0 && SWR_activity(jj,iE) == 0
                co_mat(ii,jj,iE) = -1;
                
            end
        end
    end
    %     h.CData = co_mat(:,:,iE);
    %     drawnow
    %     pause(.5)
    
    %     for iCell = length(SWR_activity(:,iE)):-1:1
    SWR_cat_data = [SWR_cat_data; this_evt];
    %     end
end


% schemaball(mean(co_mat,3)) % this doesn't work with so many cells.

% shuffle the intervals and check to see how the coactivity values change.
% Might be confounded by adjacent cells picking up each others activity.
% Could also be confounded by artifacts in the signal.

% TODO: compare the coactivity matrix for SWRs to those during increased
% REM theta and task running.




%% try using seqNMF

addpath(PARAMS.seqNMF_dir)

% SWR_cat_data =


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% old block for identifying recording epochs.
% %% Ca blocks [old: replaced with MS_extract_NLX_blocks_sandbox]
% % identify peaks in  diff(evt.t{3}) marking transitions in the camera TTLs
%
% t_idx = 3; % which event index to use.
%
% peak_threshold =  (mean(diff(nlx_evts.t{t_idx}) +0.05*std(diff(nlx_evts.t{t_idx}))));
% min_dist = 10;
% [Rec_peak, Rec_idx] = findpeaks(diff(nlx_evts.t{t_idx}), 'minpeakheight',peak_threshold, 'minpeakdistance', min_dist);
% fprintf(['\nDetected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(Rec_idx), length(Rec_idx))
%
%
% for iRec = 1:length(Rec_idx)
%     if iRec < length(Rec_idx)
%         rec_evt{iRec} = restrict(nlx_evts, nlx_evts.t{t_idx}(Rec_idx(iRec)), nlx_evts.t{t_idx}(Rec_idx(iRec+1))); % restrict the NLX evt struct to ms TTL periods
%     else
%         rec_evt{iRec} = restrict(nlx_evts, nlx_evts.t{t_idx}(Rec_idx(iRec)), nlx_evts.t{t_idx}(end)); % restrict the NLX evt file (last only)
%     end
%     all_rec_evt_len(iRec) = length(rec_evt{iRec}.t{t_idx});
% end
%
% % find the largest and use that one for now.
% [~,large_idx] = max(all_rec_evt_len);
%
%
% %% use the identified largest recording with what should be MS frame grabs
% this_csc = restrict(csc, rec_evt{large_idx}.t{t_idx}(1), rec_evt{large_idx}.t{t_idx}(end));
% fprintf('\nRestricting to section from events file. Duration: %0.0fsecs = %0.2fmins\n', (rec_evt{large_idx}.t{t_idx}(end) -rec_evt{large_idx}.t{t_idx}(1)), (rec_evt{large_idx}.t{t_idx}(end) -rec_evt{large_idx}.t{t_idx}(1))/60)
%
% %% use whole data
%
% this_csc = csc;
%
% %% initial: use a section that looks like SW sleep [use actual timestamps later but needs MS or TS files];
%
% % temp hack to test dectection
% rec.type = 'ts';
% rec.tstart = 4.413082480877258e+06; % place near the end
% rec.tend = nlx_evts.t{2}(end);
%
% this_csc = restrict(csc, rec.tstart, rec.tend);
% fprintf('\nRestricting to visually identified section.  Duration: %0.0fsecs = %0.2fmins\n',  rec.tend - rec.tstart,(rec.tend - rec.tstart)/60)


%% make a spectrogram of the average SWR

% spectrogram method using means.
% attempt to use chronux.

% %% Try Chronux?
% addpath(genpath(PARAMS.chronux_code_dir));
% disp('Chronux added to path')
% cfg_trials = [];
%
% swr_centers = IVcenters(swr_evt_out); % get the center of the swr event.
% % resize around center.
% swr_center_iv = iv([swr_centers - 0.05, swr_centers + 0.05]);
%
% % cfg_trial = []; %cfg_trial.target = this_csc.label{1}; cfg_trial.label = this_csc.label{1};
% % swr_trials = AddTSDtoIV(cfg_trial, swr_center_iv, this_csc);
%
%
% % convert data into trials
% for iEvt = length(swr_center_iv.tstart):-1:1
%
%     data_trials(iEvt,:) = this_csc.data(nearest_idx3(this_csc.data, swr_center_iv.tstart(iEvt)));
%     this_data = restrict(this_csc, swr_center_iv.tstart(iEvt), swr_center_iv.tend(iEvt));
% %     data_trials(iEvt,:) = this_data.data;
% end
%
% movingwin=[0.01 0.005]; % set the moving window dimensions
% params.Fs=this_csc.cfg.hdr{1}.SamplingFrequency; % sampling frequency
% params.fpass=[50 300]; % frequencies of interest
% params.tapers=[5 9]; % tapers
% params.trialave=1; % average over trials
% params.err=0; % no error computation
%
%
% [S1,t,f] = mtspecgramc(data_trials,movingwin,params); % compute spectrogram
%
% figure(300)
% plot_matrix(S1,t,f);
% xlabel([]); % plot spectrogram
% caxis([8 28]); colorbar;
