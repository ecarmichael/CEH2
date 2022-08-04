%% sandbox_JC_LFP_fig


% addpath(genpath('/home/williamslab/Documents/Github/CEH2'));
%
%
% inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
% addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
%
% pv1060 LTD5
% raw_dir = '/home/williamslab/Desktop/7_19_2019_PV1060_LTD5';
% decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/LTD5';
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
% ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/LTD5';
% replay_idx = [991 1001 1636 1646 2361];

% %pv1060 HATD5
% raw_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1060/11_26_2019_PV1060_HATSwitch';
% decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/HATDSwitch';
% % ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
% ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/HATDSwitch';
% replay_idx = [556 581 591 646 1206 1216 1486 1511 2291 2296 2306];

raw_dir = 'C:\Users\ecarm\Desktop';
decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1254\LTD1';
ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_13_2021_pv1254_LTD1';
replay_idx = [773 2516];

% raw_dir = 'C:\Users\ecarm\Desktop';
% decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\LTD1';
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1';
% replay_idx = [2766,3750,5847, 8171, 8336];

% raw_dir = 'C:\Users\ecarm\Desktop';
% decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1192\HATD1';
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1192\4_17_2021_PV1192_HATD1';
% replay_idx = [2584 5030 5316 ];

output_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Decoding_data\Videos\EC_museum'

re_len = 17;
core_colors = linspecer(3);


% get the files to process
f_names  = {'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\LTD1'};
replay_idx = [773, 2516];

LFP_dir = 'E:\Jisoo_Project\LFP data\Jisoo';

decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Decoding';

ft_dir = 'C:\Users\ecarm\Documents\GitHub\fieldtrip';

%% %% generate a replay event triggered spectrogram for each session
addpath(ft_dir)
ft_defaults;

for iF = 1%:length(f_names)
    
    fprintf('<strong>%s</strong>: loading ms_resize from <strong>%s</strong>\n', mfilename, f_names{iF});
    warning off; load([f_names{iF} filesep 'ms_resize.mat']); warning on;
    
    % get the session and subject IDs
    parts = strsplit(f_names{iF},  filesep);
    session = parts{end};
    subject = parts{end-1};
    date = parts{end};%(1:10);
    if strcmp(date(end), '_')
        date = date(1:end-1);
    end
    type = strsplit(parts{end}, '_');
    type = type{end};
    if strcmp(type,'HATDSwitch') && contains(subject, '1060')
        type = 'HATDSwitch';
        csc_type = 'HATSwitch';
    else
        csc_type = type;
    end
    
    % load the CSC files
    % find the right csc folder
    cd(LFP_dir)
    
    this_LFP_dir = MS_list_dir_names(cd, {subject, csc_type});
    
    cd(this_LFP_dir{1});
    
    % get the pre/post times
    EVT = LoadEvents([]);
    S_rec_idx = find(contains(EVT.label, 'Starting Recording')); % get the index for the start recoding cell
    Stop_rec_idx = find(contains(EVT.label, 'Stopping Recording')); % get the index for the start recoding cell
    
    if strcmpi(subject, 'pv1060') && strcmpi(type, 'LTD1')
        pre_S_rec_idx = 2;
        post_S_rec_idx = 3;
    elseif length(EVT.t{S_rec_idx}) ~= 2
        warning('more than two recordings detected.  Fix this later.')
        for iR = length(EVT.t{S_rec_idx}):-1:1
            rec_dur(iR) = EVT.t{Stop_rec_idx}(iR) - EVT.t{S_rec_idx}(iR);
        end
        keep_rec = find((rec_dur/60/60) > 1.5);
        if length(keep_rec) ~= 2
            error('Something is wrong with the CSC. Expected 2 recordings but found %.0f',length(keep_rec));
        end
        pre_S_rec_idx = keep_rec(1);
        post_S_rec_idx = keep_rec(2);
    else
        pre_S_rec_idx = 1;
        post_S_rec_idx = 2;
    end
    
    
    % get the channels to load from the pre-cut ms_seg data.
    cfg_load = [];
    for iC = length(ms_seg_resize.NLX_csc{1}.cfg.hdr):-1:1 % loop over channels.
        cfg_load.fc{iC} = [ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.AcqEntName '.ncs'];
        cfg_load.desired_sampling_frequency = ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.SamplingFrequency;
    end
    %     cfg_load.fc(1) = [];
    
    % hard code LFP channel
    if strcmp(subject, 'PV1043') && strcmp(type, 'LTD5')
        cfg_load.fc{1} = 'CSC6.ncs';
    end
    
    % load some data.
    for iC = length(ms_seg_resize.NLX_csc{1}.cfg.hdr):-1:1
        fprintf('<strong>%s</strong>: loading csc from <strong>%s</strong>\n', mfilename, cfg_load.fc{iC});
    end
    
    csc = MS_LoadCSC(cfg_load);
    % restrict the data to the post recording only to avoid FT issues with
    % discontinuous.
    csc = restrict(csc, EVT.t{S_rec_idx}(post_S_rec_idx), csc.tvec(end));
    
    % load the postREM binary for checking lengths of compiled post REM
    % tvec and post REM binary mat
    load([f_names{iF} filesep 'all_binary_post_REM.mat']);
    
    % get the replay events
    frame = load([decode_dir filesep subject filesep type filesep '1000shuffling.mat'], 'Final_start_frame');
    
    all_post_REM_time = [];
    for ii = 1:length(ms_seg_resize.file_names)
        if strcmpi(ms_seg_resize.hypno_label{ii}, 'REM') && strcmpi(ms_seg_resize.pre_post{ii}, 'Post')
            if length(ms_seg_resize.NLX_evt{ii}.t{end}) ~= size(ms_seg_resize.time{ii},1)
                this_time = linspace(ms_seg_resize.NLX_evt{ii}.t{end}(1), ms_seg_resize.NLX_evt{ii}.t{end}(end), size(ms_seg_resize.time{ii},1));
            else
                this_time = ms_seg_resize.NLX_evt{ii}.t{end};
            end
            all_post_REM_time = [all_post_REM_time; this_time'];
        end
    end
    
    if length(all_post_REM_time) ~= size(all_binary_post_REM, 1)
        error('compiled post REM tvec (%.0f) and saved binary (%.0f) differ in length',length(all_post_REM_time),size(all_binary_post_REM, 1) )
        
    end
    
    %% single event example
    R_idx = frame.Final_start_frame;
    win_s = 4;
    
    % R_events = findpeaks(diff(frame.Final_start_frame), 'MinPeakHeight', 15, 'MinPeakDistance', 1)
    for ii = 2:length(frame.Final_start_frame)
        if (frame.Final_start_frame(ii) - frame.Final_start_frame(ii-1)) >15
            keep_idx(ii-1) = 1;
        else
            keep_idx(ii-1) = 0;
        end
    end
    
    R_idx(~keep_idx) = [];
    
    for ii = 1:length(R_idx)
        close all
    this_event = all_post_REM_time(R_idx(ii));
    
    this_event_idx = nearest_idx3(this_event, csc.tvec);
    
    %     figure(find(contains(sessions, session))*10 + find(contains(sessions, session)))
    
%     [cfs,frq, tvec_cwt, freq_labels] = MS_plot_cw_emg(csc.tvec, csc.data(2,:), csc.data(1,:), this_event_idx, win_s);
    MS_plot_cw_emg(csc.tvec, csc.data(2,:), csc.data(1,:), this_event_idx, win_s);
   
    xlim([win_s-1 win_s+1]); % center and restrict to 1s on either side of event. Issue with labels due to cwt plot
    set(gca, 'XTickLabel', []) 
%     figure(1010)
%     subplot(2,1,1)
%     imagesc(tvec_cwt - tvec_cwt(1), frq, abs(cfs))
%     set(gca, 'YTickLabelMode', 'auto');
%     set(gca, 'YTick', freq_labels);
%     ylim([1 140])
    %         CSC_data = csc.data(1,this_event_idx - (win_s*round(1/mode(diff(csc.tvec)))):this_event_idx + (win_s*round(1/mode(diff(csc.tvec)))));
    
    %     [~,F,T,P] = spectrogram(CSC_data,rectwin(2^8), (2^8)/2, 0.5:5:80, round(1/mode(diff(csc.tvec))));
    %       ax1 = imagesc(T,F,10*log10(P));
    %                 set(ax1, 'AlphaData', ~isinf(10*log10(P)))
    %                 set(gca,'FontSize',10);
    %                 axis xy; ylabel('Frequency (Hz)');
    %                 ax = gca;
    %                 % ax.YTickLabels = {0 [], [], 60, [], [], 80}; % set
    %                 set(gca, 'tickdir', 'out');
    %                 xlim([T(1) T(end)])
    % make a raster
    
    this_ca = all_binary_post_REM(R_idx(ii) - win_s*33:R_idx(ii) + win_s*33,:)';
    this_tvec = (R_idx(ii) - win_s*33:R_idx(ii) + win_s*33)';
    this_tvec = ((this_tvec - this_tvec(1))/33) - win_s;
    
%     if exist('SA', 'var')
%         place_idx = zeros(length(this_ca),1); % allocate the index array
%         centroids = nan(size(place_idx));
%         
%         for iC = length(SA.WholePlaceCell):-1:1
%                 place_idx(SA.WholePlaceCell(iC)) = 1;
%                 centroids(SA.WholePlaceCell(iC)) = SA.PlaceFieldCentroid{SA.WholePlaceCell(iC), 3}{1}; 
%         end
%         
%         [~, cent_sort] = sort(centroids);
%         place_idx = place_idx(cent_sort);
%         
%     else
%         cent_sort = 1:size(all_binary_post,2); % just use default sort.
%         place_idx = ones(size(cent_sort));
%     end
%     
    
    
    figure(2); 
    place_idx = zeros(1,size(this_ca, 1)); 
    place_idx(SA.WholePlaceCell) = 1; 
            c_mat = [linspecer(sum(place_idx));  repmat([1 1 1], sum(~place_idx),1)]; % make colors depending on the
    MS_Ca_Raster(this_ca,this_tvec, 6,c_mat); %repmat([1,1,1], size(this_ca,1),1))
    xline(win_s, '--w', 'start', 'linewidth', 2);
    x_lim = xlim;
    xlim([-1 1])
    set(gca, 'color', 'k'); %set background color.
            colormap([linspecer(sum(place_idx));  repmat([1 1 1], 1,1)]);
    cx = colorbar;
    if exist('centroids', 'var')
        cx.TickLabels = cx.Ticks * max(centroids);
        cx.Label.String = 'place cell centroid';
    end
    
    % put the plots together
    %         axcp = copyobj(ax, fig2);
    %         set(axcp,'Position',get(ax1,'position'));
    %         delete(ax1);
    figure(1)
    figlist=get(groot,'Children');
    pause(2)
    newfig=figure(100*ii);
    tcl=tiledlayout(2,1);
    
    for jj= 1:numel(figlist)
        figure(figlist(jj));
        ax=gca;
        ax.Parent=tcl;
        ax.Layout.Tile=jj;
%         axs(jj) = gca; 
    end
    close(1);
    close(2);
    set(gcf, 'position', [662 96 758 892])
    set(gcf, 'InvertHardcopy', 'off')
    tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
    pause(2)
    if ~exist(['C:\Users\ecarm\Desktop\TEMP_LFP_figs' filesep subject filesep session], 'dir')
    mkdir(['C:\Users\ecarm\Desktop\TEMP_LFP_figs' filesep subject filesep session])
    end
    saveas(gcf, ['C:\Users\ecarm\Desktop\TEMP_LFP_figs' filesep subject filesep session filesep 'Event_' num2str(ii) '.png'])
    end
    %% event triggered spectrogram
    R_times = all_post_REM_time(frame.Final_start_frame);
    
    figure(find(contains(sessions, session))*100 + find(contains(sessions, session)))
    Triggered_Spec_FT(csc, R_times, [subject '-' type '  n = ' num2str(length(R_times))])
    hh = hline([6 12], '--r');
    hh(1).Color = [0.91 0.28 0.28];
    hh(2).Color = [0.91 0.28 0.28];
    xlim([-2 2])
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    c = colorbar;
    caxis([-8 8])
    ylabel(c, 'zscore');
    SetFigure([], gcf)
    set(gcf, 'position', [680 416 750 550])
    
    
    % event triggered LFP average
    all_lfps = [];
    win_s = [-1 1]; t_axis = win_s(1):1/csc.cfg.hdr{1}.SamplingFrequency:win_s(2);
    for iE = length(R_times):-1:1
        sta_idx = nearest_idx3([R_times(iE)+win_s(1) R_times(iE)+win_s(2)], csc.tvec);
        all_lfps(iE,:) = csc.data(1,sta_idx(1):sta_idx(2));
    end
    hold on
    plot(t_axis, nanmean(all_lfps), 'k', 'linewidth', 1.5);
    plot(t_axis, nanmean(all_lfps) + (std(all_lfps) / sqrt(size(all_lfps,1))), '--r')
    plot(t_axis, nanmean(all_lfps) - (std(all_lfps) / sqrt(size(all_lfps,1))), '--r')
    
    
end
