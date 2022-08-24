%% SCE SWR co-oc sandbox

% add vdmlad codebase and CEH2

% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2021-11-18_09-29-54_pv1252_HATD1';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2021-11-24_09-36-57_pv1252_HATDSwitch'
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2021-11-14_09-14-49_pv1252_LTD3';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-07-08_09-03-55_PV1069_LTD1';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-10-14_09-37-25_PV1069_HATD1';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-10-22_09-35-31_PV1069_HATDSwitch';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-07-12_09-24-26_PV1069_LTD5';
% LFP_data = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\LFP\2021_12_16_pv1252_MZD3_LFP';

sessions = {'E:\Jisoo_Project\LFP data\Jisoo\2021-11-18_09-29-54_pv1252_HATD1',...
    'E:\Jisoo_Project\LFP data\Jisoo\2021-11-24_09-36-57_pv1252_HATDSwitch',...
    'E:\Jisoo_Project\LFP data\Jisoo\2021-11-14_09-14-49_pv1252_LTD3',...
    'E:\Jisoo_Project\LFP data\Jisoo\2019-07-08_09-03-55_PV1069_LTD1',...
    'E:\Jisoo_Project\LFP data\Jisoo\2019-07-12_09-24-26_PV1069_LTD5',...
    'E:\Jisoo_Project\LFP data\Jisoo\2019-10-22_09-35-31_PV1069_HATDSwitch'};

ms_dir = {'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_24_2021_pv1252_HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_14_2021_pv1252_LTD3',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_12_2019_PV1069_LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch'};

Sess_pop_act_pre =[];
Sess_pop_act_post = [];
Sess_pop_shuff_pre = [];
Sess_pop_shuff_post = [];
Sess_SCE_SWR_tdiff_pre = [];
Sess_SCE_SWR_tdiff_post = [];
Sess_Shuff_SCE_SWR_tdiff_pre = [];
Sess_Shuff_SCE_SWR_tdiff_post = [];

Sess_SCE_SWR_tdiff_mean_pre = [];
Sess_SCE_SWR_tdiff_mean_post = [];

Sess_Shuff_SCE_SWR_tdiff_mean_pre = [];
Sess_Shuff_SCE_SWR_tdiff_mean_post = [];

Sess_SWR_CA_t_pre = [];
Sess_SWR_CA_t_post = [];
Sess_SCE_amp_pre = [];
Sess_SCE_amp_post = [];
Sess_Shuff_amp_pre = [];
Sess_Shuff_amp_post = [];
Sess_Shuff_SCE_SWR_t_pre = [];
Sess_Shuff_SCE_SWR_t_post = [];

plot_f = 1; % toggle plots for speed.

cutoff = .35;

for iSess = 4%:length(sessions)
    %% load some data
    
    
    LFP_data = sessions{iSess};
    
    cd(LFP_data)
    
    cfg_lfp = [];
    if contains(LFP_data, '1069')
        cfg_lfp.fc = {'CSC7.ncs'};
    else
        cfg_lfp.fc = {'CSC6.ncs'};
    end
    cfg_lfp.desired_sampling_frequency = 2000;
    csc = MS_LoadCSC(cfg_lfp);
    
    % %% some data I found in J:\Williams_Lab\Jisoo\Jisoo_Project\Inter\SWRs
    % load('J:\Williams_Lab\Jisoo\Jisoo_Project\Inter\SWRs\7_8_2019_PV1069_LTD1_activity_post.mat')
    % load('J:\Williams_Lab\Jisoo\Jisoo_Project\Inter\SWRs\7_8_2019_PV1069_LTD1_evts.mat')
    %
    %
    % length(events.SWR.tstart_idx)
    %
    % imagesc(nanmean(all_SWR_activity_post,3))
    %
    % %% look at SWRs from before
    %
    % cfg_plot = [];
    % cfg_plot.display = 'iv';
    % PlotTSDfromIV(cfg_plot, events.SWR.iv, csc)
    
    %% run the SWR detector
    cfg_swr = [];
    cfg_swr.check = 0; % plot checks.
    cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
    cfg_swr.filt.f  = [120 250]; % broad, could use 150-200?
    cfg_swr.filt.order = 4; %type filter order (fine for this f range)
    cfg_swr.filt.display_filter =0; % use this to see the fvtool
    
    
    % smoothing
    cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
    cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;
    
    % detection
    cfg_swr.threshold =2;% in sd
    cfg_swr.method = 'zscore';
    cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
    cfg_swr.merge_thr = 0.01; %merge events that are within 20ms of each other.
    %         cfg_swr.nan_idx = SWD_idx; % where are any nans, say from excluding artifacts, other events...
    
    % restrictions
    cfg_swr.max_len = [];
    cfg_swr.max_len.operation = '<';
    cfg_swr.max_len.threshold = .1;
    
    %                 cfg_swr.min_len = [];
    %                 cfg_swr.min_len.operation = '<';
    %                 cfg_swr.min_len.threshold = .2;
    cfg_swr.nCycles = 20; % number of cycles
    cfg_swr.nCycles_operation = '=<'; % number of cycles
    
    % variaence
    cfg_swr.var = [];
    cfg_swr.var.operation = '<';
    cfg_swr.var.threshold = 1;
    
    [SWR_evts, ~, SWR_amp] = MS_get_LFP_events_sandbox(cfg_swr, csc);
    
            cfg_plot.display = 'iv';
                            PlotTSDfromIV(cfg_plot, SWR_evts, csc)
    %                 pause(2)
    
    %     g = groot;
    %     if ~isempty(g.Children)
    %         figure(1)
%             saveas(gcf, 'SWR_evts_examples', 'png');
    %         close all
    %     end
    
    %% collect the events and cfgs and convert time to idx.
    parts = strsplit(cd, filesep);
    events.SWR.session = parts{end};
    events.SWR.iv = SWR_evts;
    events.SWR.cfg = cfg_swr;
    events.SWR.center_idx =  nearest_idx3(IVcenters(SWR_evts), csc.tvec);
    events.SWR.tstart_idx =  nearest_idx3(SWR_evts.tstart, csc.tvec);
    events.SWR.tend_idx =  nearest_idx3(SWR_evts.tend, csc.tvec);
    
    %         mkdir([PARAMS.inter_dir filesep 'SWRs']);
    %         save([PARAMS.inter_dir filesep 'SWRs' filesep sess_list{iSess} '_evts.mat'], 'events', '-v7.3')
    
    %% get the start and stop times for the pre and post parts of the csc.
    
    NLX_evt = LoadEvents([]);
    
    fprintf('NLX recording blocks found: %d, longest two are:  %.0d  & %.0d\n', length(NLX_evt.t{2}))
    [~, idx] = sort(NLX_evt.t{2} - NLX_evt.t{1}, 'descend');
    keep_idx = sort(idx(1:2));
    keep_idx = sort(keep_idx);
    
    pre_end = NLX_evt.t{2}(keep_idx(1));
    post_start = NLX_evt.t{1}(keep_idx(2));
    %% make a cell x time x event array.
    
    
    % move the miniscope data dir
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1')
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_24_2021_pv1252_HATDSwitch')
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_14_2021_pv1252_LTD3')
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1');
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_12_2019_PV1069_LTD5');
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch');
    % cd('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1252\2021_12_16_pv1252_MZD3');
    cd(ms_dir{iSess})
    load('ms_resize.mat');
    
    load('all_binary_pre_SW.mat');
    load('all_binary_post_SW.mat');
    
    nFrames = 33; % number of frames before and after the center of each event.
    %% get the SCEs
    
    
    % get the NLX time for Ca data
    pre_tvec = []; post_tvec = [];
    for ii = 1:length(ms_seg_resize.NLX_evt)
        if strcmp(ms_seg_resize.pre_post{ii}, 'pre') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
            pre_tvec = [pre_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
        elseif strcmp(ms_seg_resize.pre_post{ii}, 'post') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
            post_tvec = [post_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
        end
    end
    % get the indicies for recodring blocks.
    pre_seg_idx = find(diff(pre_tvec) > 1)+1;
    post_seg_idx = find(diff(post_tvec) > 1)+1;
    
    frame_n = floor(nFrames/2);
    
    if mod(frame_n, 2)==0; frame_n = frame_n+1; end
    frame_n_200 = floor(nFrames/5);
    shuff = 100;
    
    % run for pre
    data_in = all_binary_pre_SW;
    tic
    all_shuff = [];
    for iS = shuff:-1:1
        this_data = [];
        for ii = size(data_in, 2):-1:1
            this_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
        end % end cells
        
        %     all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); %malvache
        all_shuff(iS, :) = ((movsum(nansum(this_data,1), frame_n_200))./size(data_in,2))*100; % liu
    end % end shuff
    toc
    
    % find time points in real data that exceed threshold
    %  liu pop % version.
    pop_act = ((movsum(nansum(data_in,2), frame_n_200))./size(data_in,2))*100; %
    shuff_95 = prctile(all_shuff, 95, 'all');
    
    thresh = shuff_95;
    % malvache
    % thresh = shuff_mean + 3*shuff_sd;
    % shuff_mean = mean(all_shuff,'all');
    % shuff_sd = std(all_shuff, [],'all');
    % pop_act = movmean(nansum(data_in,2), frame_n_200);
    
    
    % exlude events that are too close. use findpeaks
    [peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', nFrames);
    tvec_SCE_pre =  0:(1/nFrames):(length(data_in)-1)*(1/nFrames);
    if plot_f
        figure(101);
        clf
        ax1(1) =subplot(20,1,1:7);
        MS_Ca_Raster(data_in', tvec_SCE_pre);
        % xlabel('frame number')
        ylabel('cell id')
        set(gca, 'xtick', [])
        if ~isempty(pre_seg_idx)
            vline(tvec_SCE_pre(pre_seg_idx), '--w')
        end
        ax1(2) =subplot(20,1,8:9);
        hold on
        plot(tvec_SCE_pre, pop_act);
        plot(tvec_SCE_pre(SCE_idx), pop_act(SCE_idx), 'x')
        xlim([tvec_SCE_pre(1) tvec_SCE_pre(end)])
        hline(thresh)
        ylabel('% active cells')
        set(gca, 'xtick', [])
    end
    
    % collect pre
    pop_act_pre = pop_act;
    SCE_idx_pre = SCE_idx;
    shuff_pre = all_shuff;
    thresh_pre = thresh;
    
    
    % run for post
    data_in = all_binary_post_SW;
    tic
    all_shuff = [];
    for iS = shuff:-1:1
        this_data = [];
        for ii = size(data_in, 2):-1:1
            this_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
        end % end cells
        
        %    all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); %malvache
        all_shuff(iS, :) = ((movsum(nansum(this_data,1), frame_n_200))./size(data_in,2))*100; % liu
    end % end shuff
    toc
    
    % find time points in real data that exceed threshold
    %  liu pop % version.
    pop_act = ((movsum(nansum(data_in,2), frame_n_200))./size(data_in,2))*100; %
    shuff_95 = prctile(all_shuff, 95, 'all');
    
    thresh = shuff_95;
    % malvache
    % thresh = shuff_mean + 3*shuff_sd;
    % shuff_mean = mean(all_shuff,'all');
    % shuff_sd = std(all_shuff, [],'all');
    % pop_act = movmean(nansum(data_in,2), frame_n_200);
    
    % exlude events that are too close. use findpeaks
    [peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', nFrames);
    
    tvec_SCE_post =  0:(1/nFrames):(length(data_in)-1)*(1/nFrames);
    if plot_f
        ax2(1) =subplot(20,1,11:17);
        MS_Ca_Raster(data_in', tvec_SCE_post);
        % xlabel('frame number')
        ylabel('cell id')
        set(gca, 'xtick', [])
        if ~isempty(post_seg_idx)
            vline(tvec_SCE_post(post_seg_idx), '--w')
        end
        
        ax2(2) =subplot(20,1,18:19);
        hold on
        plot(tvec_SCE_post, pop_act);
        plot(tvec_SCE_post(SCE_idx), pop_act(SCE_idx), 'x')
        xlim([tvec_SCE_post(1) tvec_SCE_post(end)])
        hline(thresh)
        ylabel('% active cells')
        set(gca, 'xtick', [])
    end
    
    % collect post
    pop_act_post = pop_act;
    SCE_idx_post = SCE_idx;
    shuff_post = all_shuff;
    thresh_post = thresh;
    
    
    
    
    pre_SCE_times = pre_tvec(SCE_idx_pre);
    post_SCE_times = post_tvec(SCE_idx_post);
        
    %%
    % initialize some arrays for later concatenation.
    all_SWR_activity_pre = [];
    all_SWR_activity_post = [];
    all_SWR_evts_pre = [];
    all_SWR_evts_post = [];
    all_SCE_activity_pre = [];
    all_SCE_activity_post = [];
    all_SCE_SWR_t_pre = [];
    all_SCE_SWR_t_post = [];
    
    all_SCE_SWR_tdiff_pre = [];
    all_SCE_SWR_tdiff_post = [];
    all_Shuff_SCE_SWR_t_diff_pre = [];
    all_Shuff_SCE_SWR_t_diff_post = [];
    
    all_Shuff_SCE_SWR_t_pre = [];
    all_Shuff_SCE_SWR_t_post = [];
    
    win_s = 3; %time in seconds
    
    for iSeg = 1:length(ms_seg_resize.NLX_csc)
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            continue
        end
        swr_r = restrict(events.SWR.iv, ms_seg_resize.NLX_csc{iSeg}.tvec(1), ms_seg_resize.NLX_csc{iSeg}.tvec(end));
        fprintf('%d %s - %s SWRs: %d\n',iSeg, ms_seg_resize.pre_post{iSeg}, ms_seg_resize.hypno_label{iSeg}, length(swr_r.tstart));
        if ~isempty(swr_r.tstart) % if there are no events here, then skip
            
            swr_centers = IVcenters(swr_r); % get the swr center.
            
            keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
            keep_idx =keep_idx(find((keep_idx ~= iSeg)));
            
            cfg_rem = []; cfg_rem.verbose = 0;
            this_ms = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
            
            this_ms = MS_de_cell(this_ms);
            
            this_ms = MS_msExtractBinary_detrendTraces(this_ms, 2);
            
            
            
            % debugging
            %         cfg_plot.display = 'tsd';
            %         cfg_plot.target = 'LFP';
            %                         ca_frames = iv(this_ms.NLX_evt.t{end}-0.001, this_ms.NLX_evt.t{end}+0.001);
            % %                         PlotTSDfromIV(cfg_plot, ca_frames, ms_seg_resize.NLX_csc{iSeg})
            %                         PlotTSDfromIV(cfg_plot, swr_r, ms_seg_resize.NLX_csc{iSeg})
            %                         PlotTSDfromIV(cfg_plot, swr_r, csc)
            
            %
            swr_idx = nearest_idx(swr_centers, this_ms.NLX_evt.t{end}); % get the SWR indicies using NLX events (TTL from frame) which correspond to the ms frame number
            
            
            count = 0; % initialize the counter.
            SWR_activity_pre = [];
            SWR_activity_post = [];
            SCE_activity_pre = [];
            SCE_activity_post = [];
            SCE_SWR_t_pre = [];
            SCE_SWR_t_post = [];
            Shuff_SCE_SWR_t_pre = [];
            Shuff_SCE_SWR_t_post = [];
            Shuff_SCE_SWR_t_diff_pre = [];
            Shuff_SCE_SWR_t_diff_post = [];
            %                     SWR_LFP_data = []; SWR_LFP_tvec = [];
            
            
            for iE = length(swr_idx):-1:1
                % check that the event is not occuring too close to the edge.
                if swr_idx(iE) <= (nFrames*win_s) || swr_idx(iE) >= length(this_ms.NLX_evt.t{end})-(nFrames*win_s)
                    continue
                else
                    count = count+1; % keep a counter for event indexing.  Avoids skipped events.
                    % get and hold the LFP values.
                    %                 this_swr = restrict(csc,this_ms.NLX_evt.t{end}(swr_idx(iE)-nFrames),this_ms.NLX_evt.t{end}(swr_idx(iE)+nFrames));
                    %                 SWR_LFP_data(count,:) = this_swr.data;
                    %                 SWR_LFP_tvec(count,:) = this_swr.tvec;
                    
                    % split into 'pre' or 'post' recording.
                    if swr_centers(iE) <= pre_end
                        for iC = this_ms.numNeurons:-1:1
                            SWR_activity_pre(iC, :, count) = this_ms.Binary(swr_idx(iE)-(nFrames*win_s): swr_idx(iE)+(nFrames*win_s), iC);
                        end
                        
                        SCE_idx = nearest_idx(swr_centers(iE), pre_tvec);
                        SCE_activity_pre(:,count) = pop_act_pre(SCE_idx-(nFrames*win_s): SCE_idx+(nFrames*win_s));
                        
                        SCE_SWR_tdiff_pre(:,count) = pre_SCE_times - swr_centers(iE);
                        
                        if sum(((pre_SCE_times - swr_centers(iE)) > -.2) & ((pre_SCE_times - swr_centers(iE)) <= 0)) >0
                            SCE_SWR_t_pre(count) = -1;
                        elseif sum(((pre_SCE_times - swr_centers(iE)) > 0) & ((pre_SCE_times - swr_centers(iE)) < .2)) >0
                            SCE_SWR_t_pre(count) = 1;
                        else
                            SCE_SWR_t_pre(count) = 0;
                        end
                        
                        % shuffle for SCE-SWR chance
                        for iShuff = 1000:-1:1
                            this_shuff = MS_randn_range(1,1,pre_SCE_times(1), pre_SCE_times(end));
                            if sum(((this_shuff - swr_centers(iE)) > -cutoff) & ((this_shuff - swr_centers(iE)) <= 0)) >0
                                Shuff_SCE_SWR_t_pre(count, iShuff) = -1;
                            elseif sum(((this_shuff - swr_centers(iE)) > 0) & ((this_shuff - swr_centers(iE)) < cutoff)) >0
                                Shuff_SCE_SWR_t_pre(count, iShuff) = 1;
                            else
                                Shuff_SCE_SWR_t_pre(count, iShuff) = 0;
                            end
                            Shuff_SCE_SWR_t_diff_post(count, iShuff) = this_shuff - swr_centers(iE);
                        end
                        
                    elseif swr_centers(iE) >= post_start
                        for iC = this_ms.numNeurons:-1:1
                            SWR_activity_post(iC, :, count) = this_ms.Binary(swr_idx(iE)-(nFrames*win_s): swr_idx(iE)+(nFrames*win_s), iC);
                        end
                        
                        SCE_idx = nearest_idx(swr_centers(iE), post_tvec);
                        SCE_activity_post(:,count) = pop_act_post(SCE_idx-(nFrames*win_s): SCE_idx+(nFrames*win_s));
                        SCE_SWR_tdiff_post(:,count) = post_SCE_times - swr_centers(iE);
                        
                        if sum(((post_SCE_times - swr_centers(iE)) > -.2) & ((post_SCE_times - swr_centers(iE)) <= 0)) >0
                            SCE_SWR_t_post(count) = -1;
                        elseif sum(((post_SCE_times - swr_centers(iE)) > 0) & ((post_SCE_times - swr_centers(iE)) < .2)) >0
                            SCE_SWR_t_post(count) = 1;
                        else
                            SCE_SWR_t_post(count) = 0;
                        end
                        
                        for iShuff = 1000:-1:1
                            this_shuff = MS_randn_range(1,1,post_SCE_times(1), post_SCE_times(end));
                            if sum(((this_shuff - swr_centers(iE)) > -cutoff) & ((this_shuff - swr_centers(iE)) <= 0)) >0
                                Shuff_SCE_SWR_t_post(count, iShuff) = -1;
                            elseif sum(((this_shuff - swr_centers(iE)) > 0) & ((this_shuff - swr_centers(iE)) < cutoff)) >0
                                Shuff_SCE_SWR_t_post(count,iShuff) = 1;
                            else
                                Shuff_SCE_SWR_t_post(count,iShuff) = 0;
                            end
                            Shuff_SCE_SWR_t_diff_post(count, iShuff) = this_shuff - swr_centers(iE);
                        end
                    end
                end
            end
        end
        % append to the master array;
        if swr_centers(end) <= pre_end
            all_SWR_activity_pre = cat(3,all_SWR_activity_pre, SWR_activity_pre);
            all_SWR_evts_pre = cat(2, all_SWR_evts_pre, swr_centers');
            all_SCE_activity_pre = cat(2,all_SCE_activity_pre, SCE_activity_pre);
            all_SCE_SWR_tdiff_pre = cat(2,all_SCE_SWR_tdiff_pre, SCE_SWR_tdiff_pre);
            all_Shuff_SCE_SWR_t_diff_pre = cat(1,all_Shuff_SCE_SWR_t_diff_pre, Shuff_SCE_SWR_t_diff_pre);
            
            all_SCE_SWR_t_pre = cat(2,all_SCE_SWR_t_pre, SCE_SWR_t_pre);
            all_Shuff_SCE_SWR_t_pre = cat(1,all_Shuff_SCE_SWR_t_pre, Shuff_SCE_SWR_t_pre);
        elseif  swr_centers(end) >= post_start
            all_SWR_activity_post = cat(3,all_SWR_activity_post, SWR_activity_post);
            all_SWR_evts_post = cat(2, all_SWR_evts_post, swr_centers');
            all_SCE_activity_post = cat(2,all_SCE_activity_post, SCE_activity_post);
            all_SCE_SWR_tdiff_post = cat(2,all_SCE_SWR_tdiff_post, SCE_SWR_tdiff_post);
            all_Shuff_SCE_SWR_t_diff_post = cat(1,all_Shuff_SCE_SWR_t_diff_post, Shuff_SCE_SWR_t_diff_post);
            
            all_SCE_SWR_t_post = cat(2,all_SCE_SWR_t_post, SCE_SWR_t_post);
            all_Shuff_SCE_SWR_t_post = cat(1,all_Shuff_SCE_SWR_t_post, Shuff_SCE_SWR_t_post);
        end
        %             all_SWR_LFP.data = cat(2,all_SWR_LFP.data, SWR_LFP_data);
        
    end % if swr_idx is empty.
    % recording segments/blocks
    
    %% Get a shuffle dist for the occurance of SWRs and SCEs within +/-200ms
    %
    % %     shuff_SCE_SWR_pre(ii) =
    %
    % rand_pre_t = MS_randn_range(10000, 1, 5, pre_end);
    %     for ii = length(rand_pre_t):-1:1
    %
    %
    %     if sum(((post_SCE_times - rand_pre_t(iE)) > -.2) & ((post_SCE_times - rand_pre_t(iE)) <= 0)) >0
    %         SCE_SWR_t_post(count) = -1;
    %     elseif sum(((post_SCE_times - rand_pre_t(iE)) > 0) & ((post_SCE_times - rand_pre_t(iE)) < .2)) >0
    %         SCE_SWR_t_post(count) = 1;
    %     else
    %         SCE_SWR_t_post(count) = 0;
    %     end
    %
    %
    % end
    
    %%  collect the SWRs from LFP
    % get the LFP snippets around the SWR
    Fs =csc.cfg.hdr{1}.SamplingFrequency; % time in seconds.
    for ii = length(all_SWR_evts_pre):-1:1
        ETAvg_pre(ii,:) = csc.data(1,nearest_idx(all_SWR_evts_pre(ii), csc.tvec)-(Fs*win_s): nearest_idx(all_SWR_evts_pre(ii), csc.tvec)+(Fs*win_s));
        SWR_Amp_Avg_pre(ii,:) = SWR_amp.data(1,nearest_idx(all_SWR_evts_pre(ii), SWR_amp.tvec)-(Fs*win_s): nearest_idx(all_SWR_evts_pre(ii), SWR_amp.tvec)+(Fs*win_s));
    end
    
    for ii = length(all_SWR_evts_post):-1:1
        ETAvg_post(ii,:) = csc.data(1,nearest_idx(all_SWR_evts_post(ii), csc.tvec)-(Fs*win_s): nearest_idx(all_SWR_evts_post(ii), csc.tvec)+(Fs*win_s));
        SWR_Amp_Avg_post(ii,:) = SWR_amp.data(1,nearest_idx(all_SWR_evts_post(ii), SWR_amp.tvec)-(Fs*win_s): nearest_idx(all_SWR_evts_post(ii), SWR_amp.tvec)+(Fs*win_s));
    end
    
    
    % get some shuffled events as well
    % fprintf('Shuffling SWR amp ')
    % tic
    % Shuff_SWR_amp_pre = [];
    % Shuff_SWR_amp_post = [];
    %
    % for iShuff = shuff:-1:1
    %     this_shuff = circshift(SWR_amp.data(1:nearest_idx3(pre_end, csc.tvec)), round(MS_randn_range(1,1,Fs*win_s, length(csc.tvec)-Fs*win_s)));
    %     Shuff_SWR_amp_pre(iShuff,:) = this_shuff(1:((Fs*win_s)*2)+1);
    %
    %     this_shuff = circshift(SWR_amp.data(nearest_idx3(post_start, csc.tvec):end), round(MS_randn_range(1,1,Fs*win_s, length(csc.tvec)-Fs*win_s)));
    %     Shuff_SWR_amp_post(iShuff,:) = this_shuff(1:((Fs*win_s)*2)+1);
    % end
    % toc
    
    
    %%  get the LFP snippets around the SCEs
    SCE_amp_pre = []; SCE_amp_post = [];
    
    for ii = length(pre_SCE_times):-1:1
        SCE_amp_pre(ii,:) = SWR_amp.data(1,nearest_idx(pre_SCE_times(ii), SWR_amp.tvec)-(Fs*win_s): nearest_idx(pre_SCE_times(ii), SWR_amp.tvec)+(Fs*win_s));
    end
    
    for ii = length(post_SCE_times):-1:1
        SCE_amp_post(ii,:) = SWR_amp.data(1,nearest_idx(post_SCE_times(ii), SWR_amp.tvec)-(Fs*win_s): nearest_idx(post_SCE_times(ii), SWR_amp.tvec)+(Fs*win_s));
    end
    
    fprintf('Shuffling SCE ripple amp ')
    tic
    Shuff_SCE_amp_pre = [];
    Shuff_SCE_amp_post = [];
    
    for iShuff = shuff:-1:1
        this_shuff = circshift(SWR_amp.data(1:nearest_idx3(pre_end, csc.tvec)), round(MS_randn_range(1,1,Fs*win_s, length(csc.tvec)-Fs*win_s)));
        Shuff_SCE_amp_pre(iShuff,:) = this_shuff(1:((Fs*win_s)*2)+1);
        
        this_shuff = circshift(SWR_amp.data(nearest_idx3(post_start, csc.tvec):end), round(MS_randn_range(1,1,Fs*win_s, length(csc.tvec)-Fs*win_s)));
        Shuff_SCE_amp_post(iShuff,:) = this_shuff(1:((Fs*win_s)*2)+1);
    end
    toc
    
    
    %% plot some SWRs, SWR-triggered average, and SWR-triggered
    % cfg_plot = [];
    % cfg_plot.display = 'iv';
    % PlotTSDfromIV(cfg_plot, Ca_SWRs, csc)
    
    tvec = -win_s:1/33:win_s;
    mean_SWR_pre = nanmean(nanmean(all_SWR_activity_pre,3),1);
    mean_SWR_post = nanmean(nanmean(all_SWR_activity_post,3),1);
    
    mean_SCE_pre = nanmean(all_SCE_activity_pre,2);
    mean_SCE_post = nanmean(all_SCE_activity_post,2);
    %% update pop x time plots
    pre_SWR_idx  = nearest_idx3(all_SWR_evts_pre, pre_tvec);
    post_SWR_idx  = nearest_idx3(all_SWR_evts_post, post_tvec);
    
    if plot_f
        figure(101);
        ax1(3) =subplot(20,1,10);
        cla
        % hold on
        line([tvec_SCE_pre(pre_SWR_idx); tvec_SCE_pre(pre_SWR_idx)], repmat([0; 1] ,1, length(pre_SWR_idx)),'color', 'r')
        set(gca, 'xtick', []);
        ylabel('SWRs')
        
        ax2(3) =subplot(20,1,20);
        cla
        line([tvec_SCE_post(post_SWR_idx); tvec_SCE_post(post_SWR_idx)], repmat([0; 1] ,1, length(post_SWR_idx)),'color', 'r')
        ylabel('SWRs')
        xlabel('Concatenated SWS time (s)')
        
        maximize
%         saveas(gcf, ['C:\Users\ecarm\Desktop\SWR_SCE_SWR_' parts{end-1} '_' parts{end} '.png'])
        
        % linkaxes(ax1, 'x')
        % linkaxes(ax2, 'x')
    end
    %%
    if plot_f
        c_ord = linspecer(4);
        figure(303)
        clf
        subplot(8,2,1)
        plot(-win_s:1/(csc.cfg.hdr{1}.SamplingFrequency):win_s, mean(ETAvg_pre),'k', 'linewidth', 2)
        ylabel('mV')
        set(gca, 'xtick', []);
        xline(0);
        title({strrep(events.SWR.session, '_', ' ' );'Pre task'})
        xlim([-1 1]);
        
        subplot(8,2,2)
        plot(-win_s:1/(csc.cfg.hdr{1}.SamplingFrequency):win_s, mean(ETAvg_post),'k', 'linewidth', 2)
        ylabel('mV')
        set(gca, 'xtick', []);
        title('Post task')
        xline(0);
        xlim([-1 1]);
        
        subplot(8,2,3:2:9)
        imagesc(-win_s:1/33:win_s, 1:size(all_SWR_activity_pre,1), mean(all_SWR_activity_pre,3))
        ylabel({'Pre' ;'Cell ID'})
        set(gca, 'xtick', []);
        xline(0, '--w', 'linewidth', 2);
        title('SWR-triggered mean Ca activity (binarized)')
        xlim([-1 1]);
        
        subplot(8,2,13:2:15)
        cla
        imagesc(-win_s:1/33:win_s,1:size(all_SCE_activity_pre,2), all_SCE_activity_pre')
        caxis([thresh_pre max(all_SCE_activity_pre, [], 'all')])
        xline(0, '--w', 'linewidth', 2);
        ylabel({ 'SWR event'})
        title('Pre SCE per SWR event')
        xlabel('Time from SWR center (s)')
        xlim([-1 1]);
        
        subplot(8,2,11)
        plot(-win_s:1/33:win_s, mean_SWR_pre*100,'color', c_ord(2,:), 'linewidth', 2)
        ylabel({'Pre' ;'mean active cells (%)'})
        ylim([min([mean_SWR_post*100, mean_SWR_pre*100]), max([mean_SWR_post*100, mean_SWR_pre*100])]);
        xline(0);
        set(gca, 'xtick', []);
        xlim([-1 1]);
        
        subplot(8,2,4:2:10)
        imagesc(-win_s:1/33:win_s, 1:size(all_SWR_activity_pre,1), mean(all_SWR_activity_post,3))
        ylabel({'Post' ;'Cell ID'})
        set(gca, 'xtick', []);
        xline(0, '--w', 'linewidth', 2);
        title('SWR-triggered mean Ca activity (binarized)')
        xlim([-1 1]);
        
        subplot(8,2,14:2:16)
        cla
        imagesc(-win_s:1/33:win_s,1:size(all_SCE_activity_post,2), all_SCE_activity_post')
        caxis([thresh_post max(all_SCE_activity_post, [], 'all')])
        xline(0, '--w', 'linewidth', 2);
        ylabel({'SWR event'})
        title('Post SCE per SWR event')
        xlabel('Time from SWR center (s)')
        xlim([-1 1]);
        
        % axis off
        % hold on
        % this_c = winter(size(all_SCE_activity_post,2));
        % for ii = 1:size(all_SCE_activity_post,2)
        %     plot(tvec, all_SCE_activity_post(:,ii)+ii*20,'color', this_c(ii,:), 'linewidth', .5);
        % end
        % ylim([min(all_SCE_activity_post(:,1)+1*10), max(all_SCE_activity_post(:,ii)+ii*10)]);
        % yline(thresh_post, '--', 'color', c_ord(4,:));
        
        subplot(8,2,12)
        plot(-win_s:1/33:win_s, mean_SWR_post*100,'color', c_ord(1,:), 'linewidth', 2)
        ylabel({'Post' ;'mean active cells (%)'})
        ylim([min([mean_SWR_post*100, mean_SWR_pre*100]), max([mean_SWR_post*100, mean_SWR_pre*100])]);
        xline(0);
        set(gca, 'xtick', []);
        xlim([-1 1])
        
        parts = strsplit(events.SWR.session, '_');
        maximize
%         saveas(gcf, ['C:\Users\ecarm\Desktop\SWR_summary_' parts{end-1} '_' parts{end} '.png'])
    end
    
    %% collect the SCE-SWR times
    all_SCE_SWR_min_pre = [];
    for ii = size(all_SCE_SWR_tdiff_pre,2):-1:1
        [~, idx] = min(abs(all_SCE_SWR_tdiff_pre(:,ii)));
        all_SCE_SWR_min_pre(ii) = all_SCE_SWR_tdiff_pre(idx,ii);
    end
    
    all_SCE_SWR_min_post = [];
    for ii = size(all_SCE_SWR_tdiff_post,2):-1:1
        [~, idx] = min(abs(all_SCE_SWR_tdiff_post(:,ii)));
        all_SCE_SWR_min_post(ii) = all_SCE_SWR_tdiff_post(idx,ii);
    end
    
    all_shuff_SCE_SWR_min_pre = [];
    for ii = size(all_Shuff_SCE_SWR_t_diff_pre,2):-1:1
        [~, idx] = min(abs(all_Shuff_SCE_SWR_t_diff_pre(:,ii)));
        all_shuff_SCE_SWR_min_pre(ii) = all_Shuff_SCE_SWR_t_diff_pre(idx,ii);
    end
    
    all_shuff_SCE_SWR_min_post = [];
    for ii = size(all_Shuff_SCE_SWR_t_diff_post,2):-1:1
        [~, idx] = min(abs(all_Shuff_SCE_SWR_t_diff_post(:,ii)));
        all_shuff_SCE_SWR_min_post(ii) = all_Shuff_SCE_SWR_t_diff_post(idx,ii);
    end
    %% follect across cells
    if exist('all_SCE_activity_pre', 'var') && ~isempty(all_SCE_activity_pre)
        Sess_pop_act_pre = [Sess_pop_act_pre; all_SCE_activity_pre'];
        Sess_pop_shuff_pre = [Sess_pop_shuff_pre; shuff_pre(:,1:size(all_SCE_activity_pre,1))];
        
        % get the time between the SCE and the closest SWR.
        Sess_SCE_SWR_tdiff_pre = [Sess_SCE_SWR_tdiff_pre, all_SCE_SWR_min_pre];
        Sess_Shuff_SCE_SWR_tdiff_pre = [Sess_Shuff_SCE_SWR_tdiff_pre, all_shuff_SCE_SWR_min_pre];
        
        Sess_SCE_SWR_tdiff_mean_pre = [Sess_SCE_SWR_tdiff_mean_pre, [nanmean(all_SCE_SWR_t_pre== 1); nanmean(all_SCE_SWR_t_pre == -1)]];
        Sess_Shuff_SCE_SWR_tdiff_mean_pre = [Sess_Shuff_SCE_SWR_tdiff_mean_pre, [nanmean(all_Shuff_SCE_SWR_t_pre== 1, 'all'); nanmean(all_Shuff_SCE_SWR_t_pre == -1, 'all')]];
        
        %keep all the time differences between SCEs and SWRs
        %         Sess_SCE_SWR_tdiff_all_pre = [Sess_SCE_SWR_tdiff_all_pre; reshape(all_SCE_SWR_tdiff_pre,[],1)];
        %         Sess_Shuff_SCE_SWR_tdiff_all_pre = [Sess_Shuff_SCE_SWR_tdiff_all_pre; reshape(all_Shuff_SCE_SWR_t_diff_pre, [],1)];
        
        Sess_SWR_CA_t_pre = [Sess_SWR_CA_t_pre; all_SCE_SWR_t_pre'];
        Sess_Shuff_SCE_SWR_t_pre = [Sess_Shuff_SCE_SWR_t_pre, all_Shuff_SCE_SWR_t_pre'];
        
    end
    if exist('all_SCE_activity_post', 'var') && ~isempty(all_SCE_activity_post)
        Sess_pop_act_post =[Sess_pop_act_post;  all_SCE_activity_post'];
        Sess_pop_shuff_post = [Sess_pop_shuff_post; shuff_post(:,1:size(all_SCE_activity_post,1))];
        
        Sess_SCE_SWR_tdiff_post = [Sess_SCE_SWR_tdiff_post, all_SCE_SWR_min_post];
        Sess_Shuff_SCE_SWR_tdiff_post = [Sess_Shuff_SCE_SWR_tdiff_post, all_shuff_SCE_SWR_min_post];
        
        
        Sess_SCE_SWR_tdiff_mean_post = [Sess_SCE_SWR_tdiff_mean_post, [nanmean(all_SCE_SWR_t_post== 1); nanmean(all_SCE_SWR_t_post == -1)]];
        Sess_Shuff_SCE_SWR_tdiff_mean_post = [Sess_Shuff_SCE_SWR_tdiff_mean_post, [nanmean(all_Shuff_SCE_SWR_t_post== 1, 'all'); nanmean(all_Shuff_SCE_SWR_t_post == -1, 'all')]];
        
        %         Sess_SCE_SWR_tdiff_all_post = [Sess_SCE_SWR_tdiff_all_post, reshape(
        
        Sess_SWR_CA_t_post = [Sess_SWR_CA_t_post; all_SCE_SWR_t_post'];
        Sess_Shuff_SCE_SWR_t_post = [Sess_Shuff_SCE_SWR_t_post all_Shuff_SCE_SWR_t_post'];
    end
        
    
    Sess_SCE_amp_pre = [Sess_SCE_amp_pre; SCE_amp_pre];
    Sess_SCE_amp_post = [Sess_SCE_amp_post; SCE_amp_post];
    
    Sess_Shuff_amp_pre = [Sess_Shuff_amp_pre; Shuff_SCE_amp_pre];
    Sess_Shuff_amp_post = [Sess_Shuff_amp_post; Shuff_SCE_amp_post];
    
    
    close all
    clearvars -except sessions ms_dir iSess win_s plot_f  Sess_* cutoff
    %_pop_act_pre Sess_pop_act_post     Sess_pop_shuff_pre Sess_pop_shuff_post    Sess_SCE_amp_pre Sess_SCE_amp_post        Sess_Shuff_amp_pre Sess_Shuff_amp_post
    
end
%% plot the population activity at the time of the SWR
c_ord = linspecer(5);
figure(106)
subplot(4,3,[1 4])
cla
hold on
% plot(-win_s:1/33:win_s, nanmean([Sess_pop_act_pre; Sess_pop_act_post]), 'b');

d_h = shadedErrorBar(-win_s:1/33:win_s, nanmean(Sess_pop_shuff_pre), nanstd(Sess_pop_shuff_pre)./sqrt(length(Sess_pop_shuff_pre)));
d_h.mainLine.Color = [.6 .6 .6];
d_h.mainLine.LineWidth = 3;
d_h.patch.FaceColor = [.6 .6 .6];
d_h.patch.FaceAlpha = .2;

d_h2 = shadedErrorBar(-win_s:1/33:win_s, nanmean(Sess_pop_act_pre), nanstd(Sess_pop_act_pre)./sqrt(length(Sess_pop_act_pre)));
d_h2.mainLine.Color = c_ord(2,:);
d_h2.mainLine.LineWidth = 3;
d_h2.patch.FaceColor = c_ord(2,:);
d_h2.patch.FaceAlpha = .2;
% control
xlabel('time from  SCE event (s)')
ylabel('Recruited cell (%)')

% plot(-win_s:1/33:win_s, nanmean([Sess_pop_shuff_pre; Sess_pop_shuff_post]), 'color', [.3 .3 .3]);
xlim([-1 1]);
xline(0);

subplot(4,3,[7 10])
cla
hold on
% control
d_h = shadedErrorBar(-win_s:1/33:win_s, nanmean(Sess_pop_shuff_post), nanstd(Sess_pop_shuff_post)./sqrt(length(Sess_pop_shuff_post)));
d_h.mainLine.Color = [.6 .6 .6];
d_h.mainLine.LineWidth = 3;
d_h.patch.FaceColor = [.6 .6 .6];
d_h.patch.FaceAlpha = .2;
xlabel('time from  SCE event (s)')
ylabel('Recruited cell (%)')

% plot(-win_s:1/33:win_s, nanmean([Sess_pop_act_pre; Sess_pop_act_post]), 'b');
d_h2 = shadedErrorBar(-win_s:1/33:win_s, nanmean(Sess_pop_act_post), nanstd(Sess_pop_act_post)./sqrt(length(Sess_pop_act_post)));
d_h2.mainLine.Color = c_ord(1,:);
d_h2.mainLine.LineWidth = 3;
d_h2.patch.FaceColor = c_ord(1,:);
d_h2.patch.FaceAlpha = .2;


% plot(-win_s:1/33:win_s, nanmean([Sess_pop_shuff_pre; Sess_pop_shuff_post]), 'color', [.3 .3 .3]);
xlim([-1 1]);
xline(0);


% get the percentage of SCE->SWR
subplot(4,3,[2])
cla
hold on
b_in = [nanmean(Sess_SCE_SWR_tdiff_mean_pre(2,:)), nanmean(Sess_SCE_SWR_tdiff_mean_pre(1,:)); nanmean(Sess_Shuff_SCE_SWR_tdiff_mean_pre(2,:)), nanmean(Sess_Shuff_SCE_SWR_tdiff_mean_pre(1,:))]'*100;
% b_h1 = bar([1 2], [nanmean(Sess_SWR_CA_t_pre == -1), nanmean(Sess_SWR_CA_t_pre == 1); nanmean(Sess_Shuff_SCE_SWR_t_pre == -1, 'all') nanmean(Sess_Shuff_SCE_SWR_t_pre == 1, 'all')]'*100);
b_h1 = bar([1 2], b_in);

b_h1(1).FaceColor = c_ord(2,:);
b_h1(2).FaceColor = [.6 .6 .6];
ylim([0 20])
set(gca, 'xtick', [1 2], 'xticklabels', {'SWR->SCE','SCE->SWR'})
ylabel('% of events')

subplot(4,3,[8])
cla
hold on
% b_in = [nanmean(Sess_SWR_CA_t_post == -1) nanmean(Sess_SWR_CA_t_post == 1); nanmean(Sess_Shuff_SCE_SWR_t_post == -1, 'all'), nanmean(Sess_Shuff_SCE_SWR_t_post == 1, 'all')]'*100;
b_in = [nanmean(Sess_SCE_SWR_tdiff_mean_post(2,:)), nanmean(Sess_SCE_SWR_tdiff_mean_post(1,:)); nanmean(Sess_Shuff_SCE_SWR_tdiff_mean_post(2,:)), nanmean(Sess_Shuff_SCE_SWR_tdiff_mean_post(1,:))]'*100;
% eb_in = [nanstd(Sess_SCE_SWR_tdiff_mean_post(2,:)), nanstd(Sess_SCE_SWR_tdiff_mean_post(1,:)); nanstd(Sess_Shuff_SCE_SWR_tdiff_mean_post(2,:)), nanstd(Sess_Shuff_SCE_SWR_tdiff_mean_post(1,:))]'*100;
b_h2 = bar([1 2], b_in);
% errorbar([1 2], b_in, eb_in)
b_h2(1).FaceColor = c_ord(1,:);
b_h2(2).FaceColor = [.6 .6 .6];
ylim([0 20])
set(gca, 'xtick', [1 2], 'xticklabels', {'SWR->SCE','SCE->SWR'})
ylabel('% of events')

subplot(4,3,[5])
cla
histogram(Sess_SCE_SWR_tdiff_pre, 25,'BinLimits', [-1 1])
% hold on
% histogram(Sess_
xlim([-1 1]);
vline([-.2 .2])
ylabel('# SCE events')
xline(-cutoff,'--k', 'SWR->SCE', 'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');
xline(cutoff,'--k', 'SCE->SWR', 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal');


subplot(4,3,[11])
cla
hold on
histogram(Sess_SCE_SWR_tdiff_post, 25,'BinLimits', [-1 1])
xlim([-1 1]);
% vline([-.2 .2])
xline(-cutoff,'--k', 'SWR->SCE', 'LabelHorizontalAlignment', 'left', 'LabelOrientation', 'horizontal');
xline(cutoff,'--k', 'SCE->SWR', 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal');

ylabel('# SCE events')
xlabel('time from SCE onset (s)')



subplot(4,3,[3 6])
cla
hold on
% control
d_h = shadedErrorBar(-win_s:1/2000:win_s, nanmean(Sess_Shuff_amp_pre), nanstd(Sess_Shuff_amp_pre)./sqrt(length(Sess_Shuff_amp_pre)));
d_h.mainLine.Color = [.6 .6 .6];
d_h.mainLine.LineWidth = 3;
d_h.patch.FaceColor = [.6 .6 .6];
d_h.patch.FaceAlpha = .2;

d_h = shadedErrorBar(-win_s:1/2000:win_s, nanmean(Sess_SCE_amp_pre), nanstd(Sess_SCE_amp_pre)./sqrt(length(Sess_SCE_amp_pre)));
d_h.mainLine.Color = c_ord(2,:);
d_h.mainLine.LineWidth = 3;
d_h.patch.FaceColor = c_ord(2,:);
d_h.patch.FaceAlpha = .2;


% plot(-win_s:1/2000:win_s, nanmean([Sess_SCE_amp_pre; Sess_SCE_amp_post]), 'r');
% plot(-win_s:1/2000:win_s, nanmean([Sess_Shuff_amp_pre; Sess_Shuff_amp_post]), 'color', [.3 .3 .3]);
xlim([-1 1]);
xline(0);
xlabel('Time from SCE (s)')
ylabel('Ripple band amplitude')

subplot(4,3,[9 12])
cla
hold on

% control
d_h = shadedErrorBar(-win_s:1/2000:win_s, nanmean(Sess_Shuff_amp_post), nanstd(Sess_Shuff_amp_post)./sqrt(length(Sess_Shuff_amp_post)));
d_h.mainLine.Color = [.6 .6 .6];
d_h.mainLine.LineWidth = 3;
d_h.patch.FaceColor = [.6 .6 .6];
d_h.patch.FaceAlpha = .2;

d_h = shadedErrorBar(-win_s:1/2000:win_s, nanmean(Sess_SCE_amp_post), nanstd(Sess_SCE_amp_post)./sqrt(length(Sess_SCE_amp_post)));
d_h.mainLine.Color = c_ord(1,:);
d_h.mainLine.LineWidth = 3;
d_h.patch.FaceColor = c_ord(1,:);
d_h.patch.FaceAlpha = .2;


% plot(-win_s:1/2000:win_s, nanmean([Sess_SCE_amp_pre; Sess_SCE_amp_post]), 'r');
% plot(-win_s:1/2000:win_s, nanmean([Sess_Shuff_amp_pre; Sess_Shuff_amp_post]), 'color', [.3 .3 .3]);
xlim([-1 1]);
xline(0);
xlabel('Time from SCE (s)')
ylabel('Ripple band amplitude')

% plot(tvec, nanmean(Sess_pop_act_post), 'b');






%% make a plot of SCE mean
% SCE_pop_act_post = [];
% 
% for ii = length(SCE_idx_pre):-1:1
% 
% 
%     SCE_pop_act_post
% 
% end
% 
% 
% figure(901)
% subplot(2,1,1)
%   imagesc(-win_s:1/33:win_s, 1:size(all_binary_pre_REM,1), mean(all_SWR_activity_pre,3))
%   ylabel({'Pre' ;'Cell ID'})
%   set(gca, 'xtick', []);
%   xline(0, '--w', 'linewidth', 2);
%   title('SWR-triggered mean Ca activity (binarized)')
%   xlim([-1 1]);
