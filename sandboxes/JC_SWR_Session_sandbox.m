%% SCE SWR co-oc sandbox

% add vdmlad codebase and CEH2

% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2021-11-18_09-29-54_pv1252_HATD1';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2021-11-24_09-36-57_pv1252_HATDSwitch'
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2021-11-14_09-14-49_pv1252_LTD3';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-07-08_09-03-55_PV1069_LTD1';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-10-14_09-37-25_PV1069_HATD1';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-10-22_09-35-31_PV1069_HATDSwitch';
% LFP_data = 'E:\Jisoo_Project\LFP data\Jisoo\2019-07-12_09-24-26_PV1069_LTD5'; 
LFP_data = '/home/williamslab/Desktop/JC_SWR/2019-07-12_09-24-26_PV1069_LTD5'; 
% LFP_data = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\LFP\2021_12_16_pv1252_MZD3_LFP';

% data_dir  =

%% load some data

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

SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc);

%         cfg_plot.display = 'tsd';
%                         PlotTSDfromIV(cfg_plot, SWR_evts, csc)
%                 pause(2)

g = groot;
if ~isempty(g.Children)
    figure(1)
    saveas(gcf, 'SWR_evts_examples', 'png');
    close all
end

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
cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_12_2019_PV1069_LTD5'); 
% cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch'); 
% cd('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1252\2021_12_16_pv1252_MZD3'); 
load('ms_resize.mat');


nFrames = 33; % number of frames before and after the center of each event.
%% get the SCEs
load('all_binary_pre_SW.mat');
load('all_binary_post_SW.mat'); 

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

    all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); 
end % end shuff
toc

shuff_mean = mean(all_shuff,'all');
shuff_sd = std(all_shuff, [],'all');

thresh = shuff_mean + 3*shuff_sd; 
% find time points in real data that exceed threshold
pop_act = movmean(nansum(data_in,2), frame_n_200);

% exlude events that are too close. use findpeaks
[peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', nFrames);

tvec_SCE =  0:(1/nFrames):(length(data_in)-1)*(1/nFrames); 
figure(101); 
ax(1) =subplot(14,1,1:5);
MS_Ca_Raster(data_in', tvec_SCE);
% xlabel('frame number')
ylabel('cell id')
set(gca, 'xtick', [])

ax(2) =subplot(14,1,6);
hold on
plot(tvec_SCE, pop_act);
plot(tvec_SCE(SCE_idx), pop_act(SCE_idx), 'x')
xlim([tvec_SCE(1) tvec_SCE(end)])
hline(thresh)
ylabel('pop activity')
set(gca, 'xtick', [])

pop_act_pre = pop_act;
SCE_idx_pre = SCE_idx;
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

    all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); 
end % end shuff
toc

shuff_mean = mean(all_shuff,'all');
shuff_sd = std(all_shuff, [],'all');

thresh = shuff_mean + 3*shuff_sd; 
% find time points in real data that exceed threshold
pop_act = movmean(nansum(data_in,2), frame_n_200);

% exlude events that are too close. use findpeaks
[peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', nFrames);

tvec_SCE =  0:(1/nFrames):(length(data_in)-1)*(1/nFrames); 
figure(101); 
ax(1) =subplot(14,1,8:12);
MS_Ca_Raster(data_in', tvec_SCE);
% xlabel('frame number')
ylabel('cell id')
set(gca, 'xtick', [])

ax(2) =subplot(13,1,13);
hold on
plot(tvec_SCE, pop_act);
plot(tvec_SCE(SCE_idx), pop_act(SCE_idx), 'x')
xlim([tvec_SCE(1) tvec_SCE(end)])
hline(thresh)
ylabel('pop activity')
set(gca, 'xtick', [])

pop_act_post = pop_act;
SCE_idx_post = SCE_idx;
thresh_post = thresh;


pre_tvec = []; post_tvec = []; 
for ii = 1:length(ms_seg_resize.NLX_evt)
    if strcmp(ms_seg_resize.pre_post{ii}, 'pre') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
        pre_tvec = [pre_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
    elseif strcmp(ms_seg_resize.pre_post{ii}, 'post') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
        post_tvec = [post_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
    end
end

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


for iSeg = 1:length(ms_seg_resize.NLX_csc)
    if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
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
        
        %                     SWR_LFP_data = []; SWR_LFP_tvec = [];
        
        for iE = length(swr_idx):-1:1
            % check that the event is not occuring too close to the edge.
            if swr_idx(iE) < nFrames || swr_idx(iE) > length(this_ms.NLX_evt.t{end})-nFrames
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
                        SWR_activity_pre(iC, :, count) = this_ms.Binary(swr_idx(iE)-nFrames: swr_idx(iE)+nFrames, iC);
                    end
                    
                    SCE_idx = nearest_idx(swr_centers(iE), pre_tvec); 
                    SCE_activity_pre(:,count) = pop_act_pre(SCE_idx-nFrames: SCE_idx+nFrames); 
                    SCE_SWR_t_pre(:,count) = pre_SCE_times - swr_centers(iE); 
                elseif swr_centers(iE) >= post_start
                    for iC = this_ms.numNeurons:-1:1
                        SWR_activity_post(iC, :, count) = this_ms.Binary(swr_idx(iE)-nFrames: swr_idx(iE)+nFrames, iC);
                    end
                    
                    SCE_idx = nearest_idx(swr_centers(iE), post_tvec);
                    SCE_activity_post(:,count) = pop_act_post(SCE_idx-nFrames: SCE_idx+nFrames);
                    SCE_SWR_t_post(:,count) = post_SCE_times - swr_centers(iE);
                end
            end
        end
        % append to the master array;
        if swr_centers(end) <= pre_end
            all_SWR_activity_pre = cat(3,all_SWR_activity_pre, SWR_activity_pre);
            all_SWR_evts_pre = cat(2, all_SWR_evts_pre, swr_centers');
            all_SCE_activity_pre = cat(2,all_SCE_activity_pre, SCE_activity_pre);
            all_SCE_SWR_t_pre = cat(2,all_SCE_SWR_t_pre, SCE_SWR_t_pre);
        elseif  swr_centers(end) >= post_start
            all_SWR_activity_post = cat(3,all_SWR_activity_post, SWR_activity_post);
            all_SWR_evts_post = cat(2, all_SWR_evts_post, swr_centers');
            all_SCE_activity_post = cat(2,all_SCE_activity_post, SCE_activity_post);
            all_SCE_SWR_t_post = cat(2,all_SCE_SWR_t_post, SCE_SWR_t_post);
        end
        %             all_SWR_LFP.data = cat(2,all_SWR_LFP.data, SWR_LFP_data);
        
    end % if swr_idx is empty.
end  % recording segments/blocks


%%  collect the SWRs from the CA recordings
for ii =  1:length(ms_seg_resize.NLX_csc)
    tstart(ii) = ms_seg_resize.NLX_csc{ii}.tvec(1);
    tend(ii) = ms_seg_resize.NLX_csc{ii}.tvec(end);
end
Ca_SWR_iv = iv(tstart, tend); 

Ca_SWRs = restrict(events.SWR.iv, Ca_SWR_iv);




Ca_SWR_cent = IVcenters(Ca_SWRs); 
ETAvg_pre = []; ETAvg_post = [];

for ii = length(all_SWR_evts_pre):-1:1
    ETAvg_pre(ii,:) = csc.data(1,nearest_idx(all_SWR_evts_pre(ii), csc.tvec)-2000*1: nearest_idx(all_SWR_evts_pre(ii), csc.tvec)+2000*1); 
end

for ii = length(all_SWR_evts_post):-1:1
    ETAvg_post(ii,:) = csc.data(1,nearest_idx(all_SWR_evts_post(ii), csc.tvec)-2000*1: nearest_idx(all_SWR_evts_post(ii), csc.tvec)+2000*1); 
end

%% plot some SWRs, SWR-triggered average, and SWR-triggered
cfg_plot = [];
cfg_plot.display = 'iv';
PlotTSDfromIV(cfg_plot, Ca_SWRs, csc)

tvec = -1:1/33:1;
mean_SWR_pre = nanmean(nanmean(all_SWR_activity_pre,3),1);
mean_SWR_post = nanmean(nanmean(all_SWR_activity_post,3),1);

mean_SCE_pre = nanmean(all_SCE_activity_pre,2);
mean_SCE_post = nanmean(all_SCE_activity_post,2);
%% update pop x time plots
% pre_SWRs = restrict(events.SWR, e

% figure(101); 
% ax(1) =subplot(14,1,1:5);
% hold on
% line(events.SWR

%%
c_ord = linspecer(4); 
figure(303)
clf
subplot(8,2,1)
plot(-1:1/2000:1, mean(ETAvg_pre),'k', 'linewidth', 2)
ylabel('mV')
set(gca, 'xtick', []); 
xline(0);
title({strrep(events.SWR.session, '_', ' ' );'Pre task'})


subplot(8,2,2)
plot(-1:1/2000:1, mean(ETAvg_post),'k', 'linewidth', 2)
ylabel('mV')
set(gca, 'xtick', []); 
title('Post task')
xline(0);

subplot(8,2,3:2:9)
imagesc(tvec, 1:size(all_SWR_activity_pre,1), mean(all_SWR_activity_pre,3))
ylabel({'Pre' ;'Cell ID'})
set(gca, 'xtick', []); 
xline(0, '--w', 'linewidth', 2);
title('SWR-triggered mean Ca activity (binarized)')

subplot(8,2,13:2:15)
cla
imagesc(tvec,1:size(all_SCE_activity_pre,2), all_SCE_activity_pre')
caxis([thresh_pre max(all_SCE_activity_pre, [], 'all')])
xline(0, '--w', 'linewidth', 2);
ylabel({ 'SWR event'})
title('Pre SCE per SWR event')

subplot(8,2,11)
plot(tvec, mean_SWR_pre,'color', c_ord(2,:), 'linewidth', 2)
ylabel({'Pre' ;'mean activity'})
ylim([min([mean_SWR_post, mean_SWR_pre]), max([mean_SWR_post, mean_SWR_pre])]);
xline(0);
set(gca, 'xtick', []); 

subplot(8,2,4:2:10)
imagesc(tvec, 1:size(all_SWR_activity_pre,1), mean(all_SWR_activity_post,3))
ylabel({'Post' ;'Cell ID'})
set(gca, 'xtick', []); 
xline(0, '--w', 'linewidth', 2);
title('SWR-triggered mean Ca activity (binarized)')

subplot(8,2,14:2:16)
cla
imagesc(tvec,1:size(all_SCE_activity_post,2), all_SCE_activity_post')
caxis([thresh_post max(all_SCE_activity_post, [], 'all')])
xline(0, '--w', 'linewidth', 2);
ylabel({'SWR event'})
title('Post SCE per SWR event')

% axis off
% hold on
% this_c = winter(size(all_SCE_activity_post,2));
% for ii = 1:size(all_SCE_activity_post,2)
%     plot(tvec, all_SCE_activity_post(:,ii)+ii*20,'color', this_c(ii,:), 'linewidth', .5);
% end
% ylim([min(all_SCE_activity_post(:,1)+1*10), max(all_SCE_activity_post(:,ii)+ii*10)]);
% yline(thresh_post, '--', 'color', c_ord(4,:));

subplot(8,2,12)
plot(tvec, mean_SWR_post,'color', c_ord(1,:), 'linewidth', 2)
ylabel({'Post' ;'mean activity'})
ylim([min([mean_SWR_post, mean_SWR_pre]), max([mean_SWR_post, mean_SWR_pre])]);
xline(0);
set(gca, 'xtick', []); 

parts = strsplit(events.SWR.session, '_'); 
saveas(gcf, ['C:\Users\ecarm\Desktop\SWR_summary_' parts{end-1} '_' parts{end} '.png'])
