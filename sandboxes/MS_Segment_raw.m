function MS_Segment_raw(cfg_in, csc_dir, ms_dir, ms_resize_dir)
%% MS_Segment_raw:
%
%
%
%    Inputs: 
%     -
%
%
%
%    Outputs: 
%     -
%
%
%
%
% EC 2020-02-18   initial version 
%
%
%%

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
cfg_csc.fc = {'CSC1.ncs','CSC7.ncs'}; % use csc files from Keys. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
cfg_csc.label = {'EMG', 'LFP'}; % custom naming for each channel.
cfg_csc.desired_sampling_frequency = 2000;
csc = MS_LoadCSC(cfg_csc); % need to comment out ExpKeys lines in LoadCSC


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

%% filtering Delta / theta

% delta filter.
cfg_filt_d = [];
cfg_filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_d.f  = [1 5];
cfg_filt_d.order = 8; %type filter order
% cfg_filt_d.display_filter = 1; % use this to see the fvtool and wait for a keypress before continuing.
delta_csc = FilterLFP(cfg_filt_d,csc);
close all


% filter into the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [6 11];
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
theta_csc = FilterLFP(cfg_filt_t, csc);

% 'wide' from Watson et al. 2016
% filter into the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [2 15];
cfg_filt_t.order = 12; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
wide_theta_csc = FilterLFP(cfg_filt_t, csc);


% get the theta_delta ratio
lfp_idx = find(strcmp(csc.label, 'LFP')); % use the LFP channel Or some other identifier in the csc.label. ATM this is to avoid using the 'EMG' channel.

td_ratio = abs(hilbert(theta_csc.data(lfp_idx, :))) ./ abs(hilbert(delta_csc.data(lfp_idx, :)));
td_ratio = filter(gausswin(4000)/sum(gausswin(4000)), 1, td_ratio); % 1d gaussian window with .5s window.

% narrow theta / wide theta (from Watson et al. 2016) [Note: that paper was
% in the neocortex.
wtd_ratio = abs(hilbert(theta_csc.data(lfp_idx, :))) ./ abs(hilbert(wide_theta_csc.data(lfp_idx, :)));
wtd_ratio = filter(gausswin(4000)/sum(gausswin(4000)), 1, wtd_ratio); % 1d gaussian window with .5s window.



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
csc.data = cat(1,csc.data ,wtd_ratio);
csc.label{end+1} = 'theta / 2-15Hz';

clear delta_csc theta_csc td_ratio wtd_ratio
fprintf('\n<strong>MS_SWR_Ca2</strong>: Delta, Theta, and Theta/Delta have been added as csc channels.\n');



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
        
%     elseif abs(length(TS{iT}.system_clock{1}) - length(evt_blocks{iT}.t{cfg_evt_blocks.t_chan})) == 1 && strcmp(cfg.TS_nlx_match, '1-sample')
%         figure(121)
%         subplot(2,1,1)
%         plot(diff(TS{iT}.system_clock{1}))
%         title(['TS timestamps. Segment: ' num2str(iT)])
%         subplot(2,1,2)
%         plot(diff(evt_blocks{iT}.t{cfg_evt_blocks.t_chan}))
%         title('NLX evt timestamps')
%         pause(1)
        
%         if length(TS{iT}.system_clock{1}) > length(evt_blocks{iT}.t{cfg_evt_blocks.t_chan})
            
        
        
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

% add in the NLX data

ms_seg = MS_append_data_sandbox(ms_seg, 'NLX_csc', res_csc, 'NLX_evt', res_evt, 'hypno_label', hypno_labels, 'time_labels', time_labels);
fprintf('\n<strong>MS_SWR_Ca2</strong>: NLX_csc appended\n');

% clear large variables from workspace for memory.
% clear ms res_csc res_evt flag




%% quick check?
check = 1; % toggle to skip check plots.
emg_chan  = find(ismember(cfg_csc.label, 'EMG')); % used to get the emg range. 

if check  ==1
    cfg_check = [];
    %     cfg_check.x_zoom = [ 0 5];
    cfg_check.Ca_type = 'RawTraces';
    %     cfg_check.Ca_type = 'BinaryTraces';
    cfg_check.chan_to_plot = ms_seg.NLX_csc{1}.label;
    cfg_check.plot_type = '3d';
    cfg_check.label = 'hypno_label';
     % get the min and max emg range for the first 5mins of the recording. used for consistency. 
    cfg_check.emg_range = [min(csc.data(emg_chan,1:(300*csc.cfg.hdr{1}.SamplingFrequency))), max(csc.data(emg_chan,1:(300*csc.cfg.hdr{1}.SamplingFrequency)))]; % get the min and max emg range for the first 10s of the recording. used for consistency. 
    
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

%% spectrogram of an episode w/
these_blocks = [3,7] %1:length(ms_seg.NLX_csc);
emg_range = [min(csc.data(emg_chan,1:(300*csc.cfg.hdr{1}.SamplingFrequency))), max(csc.data(emg_chan,1:(300*csc.cfg.hdr{1}.SamplingFrequency)))]; % get the min and max emg range for the first 10s of the recording. used for consistency. 

cut_off = 1; % do you want to use the cut_off selector?
cut_vals = NaN(2,length(ms_seg.NLX_csc)); % fill in the cutoff values.

for iBlock  = these_blocks
    win_s = 2^10; % window size for spec keep in base 2
    [~,F,T,P] = spectrogram(ms_seg.NLX_csc{iBlock}.data(2,:), rectwin(win_s), win_s/2, 0.5:.1:80, csc.cfg.hdr{1}.SamplingFrequency);
    % [~,F,T,P] = spectrogram(csc.data(2,1:4432000), rectwin(2^12), (2^12)/4, 1:.1:64, csc.cfg.hdr{1}.SamplingFrequency);
    
    
    figure(iBlock+20)
    ax_spec(1) = subplot(2,1,1);
    ax1 = imagesc(T,F,10*log10(P));
    set(ax1, 'AlphaData', ~isinf(10*log10(P)))
    set(gca,'FontSize',10);
    axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
    ax = gca;
    % ax.YTickLabels = {0 [], [], 60, [], [], 80}; % set
    set(gca, 'tickdir', 'out')
    
    % PC = pca(10*log10(P));
    % imagesc(PC(:,1))
    title([ms_seg.hypno_label{iBlock} ' block id: ' num2str(iBlock)]);
    
    ax_spec(2) = subplot(2,1,2);
    hold on
    for iChan = length(ms_seg.NLX_csc{iBlock}.label):-1:1
        this_tvec = ms_seg.NLX_csc{iBlock}.tvec-ms_seg.NLX_csc{iBlock}.tvec(1);
        if strfind(ms_seg.NLX_csc{iBlock}.label{iChan}, '/')
            hline(iChan*15, 'k')
            plot(this_tvec, ((ms_seg.NLX_csc{iBlock}.data(iChan,:)/max(ms_seg.NLX_csc{iBlock}.data(iChan,:)))*10)+iChan*15);
            text(this_tvec(1), mean(((ms_seg.NLX_csc{iBlock}.data(iChan,:)/max(ms_seg.NLX_csc{iBlock}.data(iChan,:)))*10)+iChan*15)-5, ms_seg.NLX_csc{iBlock}.label{iChan})
            
        else
            plot(this_tvec , (ms_seg.NLX_csc{iBlock}.data(iChan,:)*10000)+iChan*15);
            text(this_tvec(1), mean((ms_seg.NLX_csc{iBlock}.data(iChan,:)*10000)+iChan*15)-5, ms_seg.NLX_csc{iBlock}.label{iChan})
        end
    end
    xlim([this_tvec(1) this_tvec(end)])
    
    linkaxes(ax_spec, 'x');
%     Link = linkprop(ax_spec,{'XLim'});
%         setappdata(gcf, 'StoreTheLink', Link)
    Resize_figure
    
    if cut_off ==1 % get data points for resizing.  First point is the start, second click is the end. Don't zoom in when making the cutoffs. 
        input('Paused for inspection. Press any key to select cutoffs')
        valid_cutoff = [];
        while ~strcmp(valid_cutoff, 'y')
            [cut_x,~] = ginput(2);
            hold on
            v_ax(1) = vline(cut_x(1), 'r');
            t_ax(1) = text(cut_x(1), F(1), 'Cutoff start', 'color', 'r','FontSize',14 );
            v_ax(2) = vline(cut_x(2), 'r');
            t_ax(2) = text(cut_x(2), F(1), 'Cutoff end', 'color', 'r','FontSize',14 );
            if cut_x(1) <= this_tvec(1) && cut_x(2) >= this_tvec(end)
                valid_cutoff = input('Keep all? y/n\n', 's');
%                 keep_all = ture; 
            else
                valid_cutoff = input('Valid cutoff? y/n\n', 's');
            end
            delete(v_ax);
            delete(t_ax);
        end
%                         cut_x(1) = ms_seg.NLX_csc{iBlock}.tvec(end);
%                 cut_x(2) = ms_seg.NLX_csc{iBlock}.tvec(end);
        
        % get te 
        [~,idx]=min(abs(this_tvec-cut_x(1)));
        cut_x(1) =ms_seg.NLX_csc{iBlock}.tvec(idx);
        
        [~,idx]=min(abs(this_tvec-cut_x(2)));
        cut_x(2) =ms_seg.NLX_csc{iBlock}.tvec(idx);
        
        cut_vals(:,iBlock) = cut_x;
        
        % cut the data based on cut vals. 
        
%         ms_seg.NLX_csc{iBlock} = restrict
        
        
        
    end
    close
end


%% resize the events [WIP: has trouble resizing across ms and NLX timescales]

cfg_resize = [];
cfg_resize.tvec_to_use = 'NLX_csc'; % could be 'time', or 'NLX_csc'
cfg_resize.cutoffs = cut_vals; % should be [2 x nSegments] row 1 is start and row 2 is stop


ms_seg_resize = MS_resize_segments(cfg_resize, ms_seg); 


%% spectrogram of an episode w/
these_blocks = [3 ,7];

for iBlock  = these_blocks
    this_tvec = [];
    win_s = 2^10; % window size for spec keep in base 2
    [~,F,T,P] = spectrogram(ms_seg_resize.NLX_csc{iBlock}.data(2,:), rectwin(win_s), win_s/2, 0.5:.1:80, csc.cfg.hdr{1}.SamplingFrequency);
    % [~,F,T,P] = spectrogram(csc.data(2,1:4432000), rectwin(2^12), (2^12)/4, 1:.1:64, csc.cfg.hdr{1}.SamplingFrequency);
    
    
    figure(iBlock+200)
    ax_spec(1) = subplot(7,1,1:2);
    ax1 = imagesc(T,F,10*log10(P));
    set(ax1, 'AlphaData', ~isinf(10*log10(P)))
    set(gca,'FontSize',10, 'xtick', []);
    axis xy; xlabel('Time (s)'); ylabel('Frequency (Hz)');
    ax = gca;
    % ax.YTickLabels = {0 [], [], 60, [], [], 80}; % set
    set(gca, 'tickdir', 'out');
    xlim([T(1) T(end)])
    
    % PC = pca(10*log10(P));
    % imagesc(PC(:,1))
    title([ms_seg_resize.hypno_label{iBlock} ' block id: ' num2str(iBlock)]);
    
    ax_spec(2) = subplot(7,1,[3:6]);
    hold on
    for iChan = length(ms_seg_resize.NLX_csc{iBlock}.label):-1:1
        this_tvec = ms_seg_resize.NLX_csc{iBlock}.tvec-ms_seg_resize.NLX_csc{iBlock}.tvec(1);
        if strfind(ms_seg_resize.NLX_csc{iBlock}.label{iChan}, '/')
            hline(iChan*15, 'k')
            plot(this_tvec, ((ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)/max(ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)))*10)+iChan*15);
%             text(this_tvec(1), mean(((ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)/max(ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)))*10)+iChan*15)-5, ms_seg_resize.NLX_csc{iBlock}.label{iChan})
            yticks(iChan) = iChan*15;
            yticks_label{iChan} = ms_seg_resize.NLX_csc{iBlock}.label{iChan};
        elseif strcmp(ms_seg_resize.NLX_csc{iBlock}.label{iChan}, 'EMG') && exist('emg_range', 'var')
         ax_spec(3) = subplot(7,1,7);
         
            plot(this_tvec , ms_seg_resize.NLX_csc{iBlock}.data(iChan,:));
%             text(this_tvec(1)-2, mean(ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)), ms_seg_resize.NLX_csc{iBlock}.label{iChan})
        
            ylim(emg_range)
        else
            plot(this_tvec , (ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)*10000)+iChan*15);
%             text(this_tvec(1)-1, mean((ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)*10000)+iChan*15), ms_seg_resize.NLX_csc{iBlock}.label{iChan})
        
            yticks(iChan) = iChan*15; 
            yticks_label{iChan} = ms_seg_resize.NLX_csc{iBlock}.label{iChan};
        end
    end
    
    ax_spec(3) = subplot(7,1,7);
    y_label = get(gca, 'yticklabel');
% %     y_label = mat2cell(y_label,[1 3]);
% %     zero_idx = find(ismember(y_label, '0'));
%     y_label = (y_label(~cellfun('isempty',y_label)));
    y_label{ceil(end/2), :} = 'EMG'; 
    y_label{1} = num2str(ax_spec(3).YAxis.TickValues(1));
    y_label{end} = num2str(ax_spec(3).YAxis.TickValues(end));
    set(gca, 'yticklabel', y_label)
%     ax_spec(3).YAxis.Exponent = -3;
            xlim([this_tvec(1) this_tvec(end)])

    ax_spec(2) = subplot(7,1,[3:6]);
    yticks_label = yticks_label(~cellfun('isempty',yticks_label));
    set(gca, 'ytick', yticks,'yticklabel',yticks_label,  'xtick', [])
        xlim([this_tvec(1) this_tvec(end)])

    linkaxes(ax_spec, 'x');
    Resize_figure

end

%% clean up and export the ms_seg_resize

save([ms_resize_dir filesep 'ms_resize.mat'], 'ms_seg_resize', '-v7.3')
