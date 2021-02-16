function All_cells = Spatial_screener(cfg_in, f_info)
%% Spatial_screener_2D:
%
%
%
%    Inputs:
%    - cfg [struct]   configuration see the defaults below.
%
%    -
%
%
%
%    Outputs:
%    - All_cells [struct]  contains spatial information and stats for
%    place, speed, and acceleration.
%
%
%
%
% EC 2021-01-02   initial version
%
%
%
%% initialize

global PARAMS


cfg_def = [];
% general
cfg_def.binary_thresh =2; % number of sd for binary threshold of zscored Ca transients.
cfg_def.CI = 0.95; % confidence intervals will be equal for lower as 1-cfg.CI
cfg_def.remove_cells = []; % cells to skip for analyses. NOT USED ATM

% place
cfg_def.p_thres = 0.05; % value for pvalue cut off;
cfg_def.stability_thres = 0.5; % from van der Veldt 2020
cfg_def.nShuff = 1000;
cfg_def.p_bin_size = 3 ; % in cm
cfg_def.split_gaus_sd = 3; % sd for gaussian smoothing of place tuning for split session xcorr.

% speed
cfg_def.s_bin_size = 1.375;
cfg_def.s_bins  =  2.5:cfg_def.s_bin_size:30; % between -2cm/s^2 and 2cm/s^s with 20 bins matches van de Veldt et al. 2020
cfg_def.s_bins(cfg_def.s_bins==0) = []; %remove 0 bin.

% acceleration
cfg_def.accel_bin_size = .2;
cfg_def.accel_bins  =  -2:cfg_def.accel_bin_size:2; % between -2cm/s^2 and 2cm/s^s with 20 bins matches van de Veldt et al. 2020
cfg_def.accel_bins(cfg_def.accel_bins==0) = []; %remove 0 bin.


cfg = ProcessConfig(cfg_def, cfg_in);



%% load data

load('behav.mat')
load('ms.mat')
ms = MS_msExtractBinary_detrendTraces(ms, cfg.binary_thresh);

% cfg_rm.remove_idx = [17 26];
% ms = MS_Remove_trace(cfg_rm, ms);

figure(101)
% subplot(4,4,[5,6,7, 9,10,11, 13,14,15])
cfg_plot = [];
cfg_plot.view =[0 75];
cfg_plot.plot_type = '2d';
% cfg_plot.colors = parula(size(ms.Binary,2));
MS_plot_ca(cfg_plot, ms)
xlabel('time (s)')
title([f_info.subject ' ' f_info.date ' ' f_info.task])

%% align behaviour and Ca

% check if the behav needs to be interpolated.
if behav.time(end) ~= ms.time(end) || length(behav.time) ~= length(ms.time)
    fprintf('<strong> %s </strong>: behaviour and Ca are not the same length or end time.  attempting alignment \n', mfilename);
    
    behav_aligned = MS_align_data(behav, ms);
else 
    behav_aligned = behav; 
end

%smooth speed
behav_aligned.speed = smooth(behav_aligned.speed, 3*mode(diff(ms.time)));

movement_idx = behav_aligned.speed >cfg.s_bins(1); % get times when the animal was moving.
accel_movement_idx = behav_aligned.speed(1:end-1) >cfg.s_bins(1); % same as above but correct for diff used in acceleration calculation.

% get the acceleration
behav_aligned.accel = diff(smooth(behav_aligned.speed, 3*mode(diff(ms.time))));



%% plot basics for each cell

for iC = 1:size(ms.Binary,2) % loop through cells
    fprintf('\nProcessing cell %d...', iC);
    %% activity measures
    
    % number of transients
    ca_evts = MS_get_binary_events(ms.Binary(:,iC));
    
    % activity index
    diff_act = diff(ms.Binary(:,iC)); % get the change in activity states
    
    activity_idx = sum(diff_act == 1)/length(diff_act); % get the p (inactive -> active)
    
    %         % bursting index p(active -> active)  % TODO make this A) work B)
    %         more efficient
    %         for ii = length(ms.Binary(:,iC)):-1:2
    %            if ms.Binary(ii-1,iC) == 1  && ms.Binary(ii,iC) == 1
    %             burst_tran(ii) = 1;
    %            else
    %                burst_tran(ii) = 0;
    %            end
    %         end
    %
    %
    %         for ii = length(x):-1:2
    %            if x(ii-1) == 1  && x(ii) == 1
    %             burst_tran(ii) = 1;
    %            else
    %                burst_tran(ii) = 0;
    %            end
    %         end
    
    %% get the place information and stats
    
    X_bins = min(behav_aligned.position(:,1)):cfg.p_bin_size:max(behav_aligned.position(:,1));
    X_bin_centers = X_bins +  cfg.p_bin_size/2;
    X_bin_centers = X_bin_centers(1:end-1);
    % same for Y bins
    Y_bins = min(behav_aligned.position(:,2)):cfg.p_bin_size:max(behav_aligned.position(:,2));
    Y_bin_centers = Y_bins +  cfg.p_bin_size/2;
    Y_bin_centers = Y_bin_centers(1:end-1);
    
    % get place information
    [Place_MI, Place_posterior, Place_occupancy, ~, Place_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx,iC),behav_aligned.position(movement_idx,1), X_bins);
    
    % get shuffle data
    Place_shuff_tuning_curve = MS_split_shuff(ms.Binary(:,iC), behav_aligned.position(:,1),movement_idx, cfg.nShuff, X_bins);
    
    % get stats
    Place_stats= MS_boot_shuff(ms.Binary(:,iC), behav_aligned.position,movement_idx, cfg.nShuff, X_bins, Y_bins);
    
    %get sig tuning
    pval = sum(Place_shuff_tuning_curve > Place_tuning_curve,3)/cfg.nShuff;
    Place_Sig_map = Place_tuning_curve;
    Place_Sig_map(pval > cfg.p_thres) = 0;
    
    
    % split session stats and info.
    split_1 = zeros(size(ms.Binary(:,1)));
    
    %     % hard split based on time.
    %     split_1(1:ceil(length(ms.Binary(:,1))/2)) = 1;
    
    % split based on number of transients;
    split_evt_idx = ca_evts(ceil(length(ca_evts)/2),1); % get the start of the middle event.
    split_1(1:split_evt_idx) = 1;
    
    % make keep indices for split halves.
    split_1 = logical(split_1);
    split_2 = logical(~split_1);
    
    % split 1
    [Place_S1_MI, ~, Place_S1_occupancy, ~, Place_S1_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx & split_1,iC),behav_aligned.position(movement_idx & split_1,1), X_bins);
    Place_S1_shuff_tuning_curve = MS_split_shuff(ms.Binary(split_1,iC), behav_aligned.position(split_1,1),movement_idx, cfg.nShuff, X_bins);
    pval_S1 = sum(Place_S1_shuff_tuning_curve > Place_S1_tuning_curve,3)/cfg.nShuff;
    
    Place_S1_Sig_map = Place_S1_tuning_curve;
    Place_S1_Sig_map(pval_S1 > cfg.p_thres) = 0;
    
    % split 2
    [Place_S2_MI, ~, Place_S2_occupancy, ~, Place_S2_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx & split_2,iC),behav_aligned.position(movement_idx & split_2,1), X_bins);
    Place_S2_shuff_tuning_curve = MS_split_shuff(ms.Binary(split_2,iC), behav_aligned.position(split_2,1),movement_idx, cfg.nShuff, X_bins);
    pval_S2 = sum(Place_S2_shuff_tuning_curve > Place_S2_tuning_curve,3)/cfg.nShuff;
    
    Place_S2_Sig_map = Place_S2_tuning_curve;
    Place_S2_Sig_map(pval_S2 > cfg.p_thres) = 0;
    
    
    % smooth with guassian
    S1_tuning_curve_smooth = smooth(Place_S1_tuning_curve, 2);
    S2_tuning_curve_smooth = smooth(Place_S2_tuning_curve, 2);
    Place_Stability_corr = corr(S1_tuning_curve_smooth, S2_tuning_curve_smooth);
    
    % get the 95% CI for split shuff corr for comparison to real data
    for iShuff = cfg.nShuff:-1:1
        S1_shuff_tuning_curve_smooth = smooth(Place_S1_shuff_tuning_curve(:,iShuff), 2);
        S2_shuff_tuning_curve_smooth = smooth(Place_S2_shuff_tuning_curve(:,iShuff), 2);
        Place_Shuff_Stability_corr(iShuff) = corr(S1_shuff_tuning_curve_smooth, S2_shuff_tuning_curve_smooth);
    end
    
    % GE method
    sorted_Place_shuff_corr = sort(Place_Shuff_Stability_corr);
    CI_idx_loc = cfg.CI*cfg.nShuff/2;
    median_idx = round(cfg.nShuff/2);
    upper_CI95_idx = median_idx+CI_idx_loc;
    lower_CI95_idx = median_idx-CI_idx_loc;
    
    % This will make sure that upper and lower bounds are withing the actual bootstrap data
    upper_CI95_idx(upper_CI95_idx > cfg.nShuff) = cfg.nShuff;
    upper_CI95_idx(upper_CI95_idx < 1) = 1;
    lower_CI95_idx(lower_CI95_idx > cfg.nShuff) = cfg.nShuff;
    lower_CI95_idx(lower_CI95_idx < 1) = 1;
    
    Place_split_stats.upper_CI95 = sorted_Place_shuff_corr(upper_CI95_idx); % get upp CI
    Place_split_stats.lower_CI95 = sorted_Place_shuff_corr(lower_CI95_idx); % get lower CI
    
    Place_split_stats.percentile_95 = prctile(Place_Shuff_Stability_corr, cfg.CI*100); % add the percentile for good measure
    
    
    % at least one significant spatial bin
    if ~isempty(Place_Sig_map(Place_Sig_map ~=0)) % any sig spatial points
        Place_map_sig = 1; else; Place_map_sig = 0;
    end
    %
    %     % check for MI sig
    %     if mean(Place_stats.actual_MI) >   mean(Place_stats.upper_CI95) % passes MI 95CI of shuffle data
    %         Place_MI_sig = 1; else;  Place_MI_sig = 0;
    %     end
    %
    % check for split corr sig
    if Place_Stability_corr >   Place_split_stats.upper_CI95 % passes split stability 95CI of shuffle data
        Place_split_sig = 1; else; Place_split_sig = 0;
    end
    
    %% speed information
    % make bins
    Speed_bin_centers = cfg.s_bins + cfg.s_bin_size/2;
    Speed_bin_centers = Speed_bin_centers(1:end-1);
    
    %get the spatial info
    [Speed_MI, ~, Speed_occupancy,Speed_p_active, Speed_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx,iC), behav_aligned.speed(movement_idx), cfg.s_bins);
    
    % get shuffle data
    Speed_shuff_tuning_curve = MS_split_shuff(ms.Binary(movement_idx,iC), behav_aligned.speed(movement_idx),movement_idx, cfg.nShuff, cfg.s_bins);
    
    %get stats
    Speed_stats= MS_boot_shuff(ms.Binary(movement_idx(1:end-1),iC), behav_aligned.speed,movement_idx(1:end-1), cfg.nShuff, cfg.s_bins);
    
    %get sig tuning
    pval = sum(Speed_shuff_tuning_curve > Speed_tuning_curve,2)/cfg.nShuff;
    
    Speed_Sig_map = Speed_tuning_curve;
    Speed_Sig_map(pval > cfg.p_thres) = 0;
    
    %     % split session stats and info.
    %     [Speed_S1_MI, ~, Speed_S1_occupancy,Speed_S1_p_active, Speed_S1_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx & split_1,iC),ms.time(movement_idx & split_1), behav_aligned.speed(movement_idx & split_1), cfg.s_bins);
    %     Speed_S1_shuff_tuning_curve = MS_split_shuff(ms.Binary(movement_idx & split_1,iC), behav_aligned.speed(movement_idx & split_1),movement_idx, cfg.nShuff,cfg.s_bins);
    %     Speed_S1_Sig_TC = sum(Speed_S1_shuff_tuning_curve > Speed_S1_tuning_curve,3)/cfg.nShuff;
    %
    %     Speed_S1_Sig_map = Speed_S1_tuning_curve;
    %     Speed_S1_Sig_map(Speed_S1_Sig_TC < cfg.p_thres) = 0;
    %
    %
    %     [Speed_S2_MI, ~, Speed_S2_occupancy,Speed_S2_p_active, Speed_S2_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx & split_2,iC),ms.time(movement_idx & split_2), behav_aligned.speed(movement_idx & split_2), cfg.s_bins);
    %     Speed_S2_shuff_tuning_curve = MS_split_shuff(ms.Binary(movement_idx & split_2,iC), behav_aligned.speed(movement_idx & split_2),movement_idx, cfg.nShuff,cfg.s_bins);
    %     Speed_S2_Sig_TC = sum(Speed_S2_shuff_tuning_curve > Speed_S2_tuning_curve,3)/cfg.nShuff;
    %
    %     Speed_S2_Sig_map = Speed_S2_tuning_curve;
    %     Speed_S2_Sig_map(Speed_S2_Sig_TC < cfg.p_thres) = 0;
    %
    %     % smooth with guassian
    %     Speed_S1_tuning_curve_smooth = imgaussfilt(Speed_S1_tuning_curve, cfg.split_gaus_sd);
    %     Speed_S2_tuning_curve_smooth = imgaussfilt(Speed_S2_tuning_curve, cfg.split_gaus_sd);
    %
    %     Speed_Stability_corr = corr2(Speed_S1_tuning_curve_smooth, Speed_S2_tuning_curve_smooth);
    %
    %     % get the 95% CI for split shuff corr for comparison to real data
    %     for iShuff = cfg.nShuff:-1:1
    %         S1_shuff_tuning_curve_smooth = imgaussfilt(Speed_S1_tuning_curve_smooth(:,iShuff), 2);
    %         S2_shuff_tuning_curve_smooth = imgaussfilt(Speed_S2_tuning_curve_smooth(:,iShuff), 2);
    %         Speed_Shuff_Stability_corr(iShuff) = corr2(S1_shuff_tuning_curve_smooth, S2_shuff_tuning_curve_smooth);
    %     end
    %
    %     % GE method
    %     sorted_Speed_shuff_corr = sort(Speed_Shuff_Stability_corr,2);
    %     CI_idx_loc = cfg.CI*cfg.nShuff/2;
    %     median_idx = round(cfg.nShuff/2);
    %     upper_CI95_idx = median_idx+CI_idx_loc;
    %     lower_CI95_idx = median_idx-CI_idx_loc;
    %
    %     % This will make sure that upper and lower bounds are withing the actual bootstrap data
    %     upper_CI95_idx(upper_CI95_idx > cfg.nShuff) = cfg.nShuff;
    %     upper_CI95_idx(upper_CI95_idx < 1) = 1;
    %     lower_CI95_idx(lower_CI95_idx > cfg.nShuff) = cfg.nShuff;
    %     lower_CI95_idx(lower_CI95_idx < 1) = 1;
    
    %     Speed_split_stats.upper_CI95 = sorted_Speed_shuff_corr(:,upper_CI95_idx); % get upp CI
    %     Speed_split_stats.lower_CI95 = sorted_Speed_shuff_corr(:,lower_CI95_idx); % get lower CI
    %
    %     Speed_split_stats.percentile_95 = prctile(Speed_Shuff_Stability_corr, cfg.CI*100); % add the percentile for good measure
    %
    %
    % at least one significant spatial bin
    if ~isempty(Speed_Sig_map(Speed_Sig_map ~=0)) % any sig spatial points
        Speed_map_sig = 1; else; Speed_map_sig = 0;
    end
    
    % check for MI sig
    if Speed_MI >   Speed_stats.upper_CI95 % passes MI 95CI of shuffle data
        Speed_MI_sig = 1; else;  Speed_MI_sig = 0;
    end
    
    %     % check for split corr sig
    %     if Speed_Stability_corr >   Speed_split_stats.upper_CI95 % passes split stability 95CI of shuffle data
    %         Speed_split_sig = 1; else; Speed_split_sig = 0;
    %     end
    %
    
    
    %% acceleration information
    
    Acc_bin_centers = cfg.accel_bins + cfg.accel_bin_size/2;
    Acc_bin_centers = Acc_bin_centers(1:end-1);
    
    % for debugging.
    %     accel_hist = hist(behav_aligned.accel(movement_idx(1:end-1)), A_bin_centers);
    %     accel_active_hist = hist(behav_aligned.accel(ms.Binary(1:end-1,iC) & movement_idx(1:end-1)), A_bin_centers);
    %     yyaxis left
    %     bar(A_bin_centers, accel_hist./mode(diff(behav_aligned.time)), 'facecolor', PARAMS.D_grey, 'edgecolor', PARAMS.D_grey)
    %     ylabel('time in acceleration (s)')
    %
    %     %         plot(S_bin_centers, speed_active_hist./mode(diff(behav_aligned.time)), 'color', PARAMS.red)
    %     yyaxis right
    %     plot(A_bin_centers, accel_active_hist./accel_hist, 'color', PARAMS.red)
    %     %     legend({'occupancy', 'acitve'});
    %     xlabel('acceleration (cm/s^2)')
    %     ylabel('p active')
    
    % get acceleration information
    [Acc_MI, ~, Acc_occupancy,Acc_p_active, Acc_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx(1:end-1),iC),behav_aligned.accel(movement_idx(1:end-1)), cfg.accel_bins);
    
    % get shuffle data
    Acc_shuff_tuning_curve = MS_split_shuff(ms.Binary(movement_idx(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1)),movement_idx(1:end-1), cfg.nShuff, cfg.accel_bins);
    % get stats
    Acc_stats= MS_boot_shuff(ms.Binary(movement_idx(1:end-1),iC), behav_aligned.accel,movement_idx(1:end-1), cfg.nShuff, cfg.accel_bins);
    %get sig tuning
    Acc_Sig_TC = sum(Acc_shuff_tuning_curve > Acc_tuning_curve,2)/cfg.nShuff;
    
    Acc_Sig_map = Acc_tuning_curve;
    Acc_Sig_map(Acc_Sig_TC < cfg.p_thres) = 0;
    
    % split session stats and info.
    [Acc_S1_MI, ~, Acc_S1_occupancy,Acc_S1_p_active, Acc_S1_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx(1:end-1) & split_1(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1) & split_1(1:end-1)), cfg.accel_bins);
    Acc_S1_shuff_tuning_curve = MS_split_shuff(ms.Binary(movement_idx(1:end-1) & split_1(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1) & split_1(1:end-1)),movement_idx, cfg.nShuff,cfg.accel_bins);
    Acc_S1_Sig_TC = sum(Acc_S1_shuff_tuning_curve > Acc_S1_tuning_curve,3)/cfg.nShuff;
    
    Acc_S1_Sig_map = Acc_S1_tuning_curve;
    Acc_S1_Sig_map(Acc_S1_Sig_TC < cfg.p_thres) = 0;
    
    
    [Acc_S2_MI, ~, Acc_S2_occupancy,Acc_S2_p_active, Acc_S2_tuning_curve] = MS_get_spatial_information(ms.Binary(movement_idx(1:end-1) & split_2(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1) & split_2(1:end-1)), cfg.accel_bins);
    Acc_S2_shuff_tuning_curve = MS_split_shuff(ms.Binary(movement_idx(1:end-1) & split_2(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1) & split_2(1:end-1)),movement_idx, cfg.nShuff,cfg.accel_bins);
    Acc_S2_Sig_TC = sum(Acc_S2_shuff_tuning_curve > Acc_S2_tuning_curve,3)/cfg.nShuff;
    
    Acc_S2_Sig_map = Acc_S2_tuning_curve;
    Acc_S2_Sig_map(Acc_S2_Sig_TC < cfg.p_thres) = 0;
    
    
    % smooth with guassian
    Acc_S1_tuning_curve_smooth = imgaussfilt(Acc_S1_tuning_curve, cfg.split_gaus_sd);
    Acc_S2_tuning_curve_smooth = imgaussfilt(Acc_S2_tuning_curve, cfg.split_gaus_sd);
    
    Acc_Stability_corr = corr2(Acc_S1_tuning_curve_smooth, Acc_S2_tuning_curve_smooth);
    
    
    % compute split corr for shuffle and use for comparing to real
    % split corr.  If >95% CI then classify as a place cell.
    
    % if the number of transients is < 3 then this will not be
    % realistic.
    
    % if statement.
    
    %% collect information from each cell.
    All_cells.fname{iC} = f_info;

    % place info
    All_cells.place.MI(iC) = Place_MI;
    All_cells.place.Sig_map(:,iC) = Place_Sig_map;
    All_cells.place.occupanyc(:,iC) = Place_occupancy;
    All_cells.place.stats{iC} =  Place_stats;
    All_cells.place.n_evts(iC) = length(ca_evts);
    All_cells.place.ca_evts{iC} = ca_evts;
    All_cells.place.map_sig(iC) = Place_map_sig;
    
    
    All_cells.place.split.Stability_corr(iC) = Place_Stability_corr;
    All_cells.place.split.split_sig(iC) = Place_split_sig;
    All_cells.place.split.stats = Place_split_stats;
    
    All_cells.place.split.S1_MI(iC) = Place_S1_MI;
    All_cells.place.split.S1_Sig_map(:,iC) = Place_S1_Sig_map;
    All_cells.place.split.S1_occupanyc(:,iC) = Place_S1_occupancy;
    
    All_cells.place.split.S2_MI(iC) = Place_S2_MI;
    All_cells.place.split.S2_Sig_map(:,iC) = Place_S2_Sig_map;
    All_cells.place.split.S2_occupanyc(:,iC) = Place_S2_occupancy;
    
    
    
    % speed
    All_cells.speed.MI(iC) = Speed_MI;
    All_cells.speed.Sig_map(:,iC) = Speed_Sig_map;
    All_cells.speed.occupanyc(:,iC) = Speed_occupancy;
    All_cells.speed.stats{iC} = Speed_stats;
    
    %     All_cells.speed.split.S1_MI(iC) = Speed_S1_MI;
    %     All_cells.speed.split.S1_Sig_map(:,:,iC) = Speed_S1_Sig_map;
    %     All_cells.speed.split.S1_occupanyc(:,:,iC) =Speed_occupancy;
    %
    %     All_cells.speed.split.S2_MI(iC) = Speed_S2_MI;
    %     All_cells.speed.split.S2_Sig_map(:,:,iC) = Speed_S2_Sig_map;
    %     All_cells.speed.split.S2_occupanyc(:,:,iC) = Speed_S2_occupancy;
    %     All_cells.speed.split.Stability_corr(iC) = Speed_Stability_corr;
    
    
    % acceleration
    All_cells.accel.MI(iC) = Acc_MI;
    All_cells.accel.Sig_map(:,iC) = Acc_Sig_map;
    All_cells.accel.occupanyc(:,iC) = Acc_occupancy;
    All_cells.accel.stats{iC} = Acc_stats;
    
    %     All_cells.accel.split.S1_MI(iC) = Acc_S1_MI;
    %     All_cells.accel.split.S1_Sig_map(:,:,iC) = Acc_S1_Sig_map;
    %     All_cells.accel.split.S1_occupanyc(:,:,iC) =Acc_occupancy;
    %
    %     All_cells.accel.split.S2_MI(iC) = Acc_S2_MI;
    %     All_cells.accel.split.S2_Sig_map(:,:,iC) = Acc_S2_Sig_map;
    %     All_cells.accel.split.S2_occupanyc(:,:,iC) = Acc_S2_occupancy;
    %     All_cells.accel.split.Stability_corr(iC) = Acc_Stability_corr;
    
    
    %% figure 2 place information
    figure(iC)
    
    subplot(6,4,[1 5])
    axis off
    text(0,1, ['Cell # ' num2str(iC)]);
    text(0,.66,'Whole session')
    text(0,.33,['MI: ' num2str(Place_MI,2)]);
    text(0,0, ['Num trans: ' num2str(length(ca_evts))]);
    colormap(gca, 'cool')
    colorbar('location', 'southoutside', 'ticks', [0, 1], 'ticklabels', {'1^s^t', 'last'})
    
    
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,2)
        axis off
        title('Activity')
        subplot(6,4,6)
    else
        subplot(6,4,[2 6])
    end
    t_binary = ms.Binary(:,iC) & movement_idx;
    hold on
    plot(behav_aligned.time, behav_aligned.position(:,1), 'color', PARAMS.L_grey)
    set(gca, 'xtick', [], 'ytick', []);
    % put dots on positions when the cell was active.
    MS_color_plot(behav_aligned.time(t_binary), behav_aligned.position(t_binary,1), '.', cool(length(behav_aligned.position(t_binary,1))))
    axis off

    
    % overall occupancy map
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,3)
        title('Occupancy')
        axis off
        subplot(6,4,7);
    else
        subplot(6,4,[3 7]);
        title('Occupancy')
    end
   imagesc('XData', X_bin_centers, 'CData', Place_occupancy');
    set(gca, 'xtick', [], 'ytick', []);
    axis off
    
    % overall tuning map
    % overall occupancy map
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,4)
        title(['Sig at p < ' num2str(cfg.p_thres)])
        axis off
        subplot(6,4,8);
    else
        subplot(6,4,[4 8]);
        title(['Sig at p < ' num2str(cfg.p_thres)])
    end
   imagesc('XData', X_bin_centers, 'CData',Place_Sig_map');
    hold on
    set(gca, 'xtick', [], 'ytick', []);
    axis off
   
    
    % first half
    
    subplot(6,4,[9 13])
    axis off
    text(0,.8,'1^s^t half split')
    text(0,.6,['MI: ' num2str(Place_S1_MI,2)]);
        
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,14)
    else
        subplot(6,4,[10 14])
    end
    t_binary = ms.Binary(:,iC) & movement_idx & split_1;
    hold on
    plot(behav_aligned.time(movement_idx & split_1), behav_aligned.position(movement_idx & split_1,1), 'color', PARAMS.L_grey)

    set(gca, 'xtick', [], 'ytick', []);
    % put dots on positions when the cell was active.
    plot(behav_aligned.time(t_binary,1), behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red)
    axis off

    % overall occupancy map
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,15);
    else
        subplot(6,4,[11 15]);
    end
    imagesc('XData', X_bin_centers, 'CData',  Place_S1_occupancy');
    set(gca, 'xtick', [], 'ytick', []);
    axis off

    % 1st half tuning map
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,16);
    else
        subplot(6,4,[12 16]);
    end
    imagesc('XData', X_bin_centers, 'CData', Place_S1_tuning_curve');
    set(gca, 'xtick', [], 'ytick', []);
    axis off

    
    % second half
    
    subplot(6,4,17)
    axis off
    text(0,.8,'2^n^d half split')
    text(0,.6,['MI: ' num2str(Place_S2_MI,2)]);
    text(0,.4,['split xcorr: ' num2str(Place_Stability_corr,2)]);
    
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,22)
    else
        subplot(6,4,[18 22])
    end
    t_binary = ms.Binary(:,iC) & movement_idx & split_2;
    hold on
    plot(behav_aligned.time(movement_idx & split_2), behav_aligned.position(movement_idx & split_2,1), 'color', PARAMS.L_grey)

    set(gca, 'xtick', [], 'ytick', []);
    % put dots on positions when the cell was active.
    plot(behav_aligned.time(t_binary), behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red)
    axis off

    % 2nd half occupancy
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,23);
    else
        subplot(6,4,[19 23]);
    end
    imagesc('XData', X_bin_centers, 'CData', Place_S2_occupancy');
    set(gca, 'xtick', [], 'ytick', []);
    axis off
    
    % overall tuning map
    if contains(All_cells.fname{iC}.task, 'rec') || contains(All_cells.fname{iC}.task, 'LT')
        subplot(6,4,24);
    else
        subplot(6,4,[20 24]);
    end
    hold on
    imagesc('XData', X_bin_centers, 'CData', Place_S2_tuning_curve');
    set(gca, 'xtick', [], 'ytick', []);
    axis off
    x_lim = xlim;
    y_lim = ylim;
    plot([x_lim(2)-8, x_lim(2)],[y_lim(1)-(.1*y_lim(1)), y_lim(1)-(.1*y_lim(1))], 'k', 'linewidth', 1)
    text(x_lim(2)-8, y_lim(1)-(.5*y_lim(1)), '10cm', 'fontsize', 6)
%     
    saveas(gcf, [PARAMS.inter_dir 'Place_figs' filesep f_info.subject '_' f_info.date '_' f_info.task '_Cell_' num2str(iC)], 'png')
    saveas(gcf, [PARAMS.inter_dir  'Place_figs' filesep f_info.subject '_' f_info.date '_' f_info.task '_Cell_' num2str(iC)], 'fig')
%     
    %% plot everything
    %     if ishandle(300)
    %         close(300)
    %     end
    figure(100+iC)
    
    M = 4; % rows
    N = 5; % columns
    %   fig{iC} = figure('Visible', 'off'); % hack to stop figures from taking
    %   over the screen.  Good for batch processing in the background.
    
    %get binary 'event times' to be plotted as dots
    t_binary = ms.Binary(:,iC) & movement_idx;
%     accel_t_binary = find(ms.Binary(1:end-1,iC)==1);
    accel_t_binary = t_binary(1:end-1); 
    
    
    %%% title information
    subplot(M, N, 4:5) % title information. Top right corner,
    ylim([0 10])
    text(0,9,['Cell id: ' num2str(iC)], 'fontweight', 'bold')
    text(0,7,['Subject: ' f_info.subject]);
    text(0,5,['Session: ' f_info.task]);
    text(0,3,['Date: ' f_info.date]);
    text(0,1,['Binary thresh: ' num2str(ms.Binary_threshold) 'sd'  '    Num transients: ' num2str(length(ca_evts))])
    
    axis off
    
    
    %%% raw trace
    subplot(M, N, 1:2)
    plot(ms.time/1000, ms.RawTraces(:,iC), 'color', PARAMS.blue)
    xlim([ms.time(1)/1000 ms.time(end)/1000]);
    xlabel('time(s)');
    ylabel('dF/F');
    hline(mean(ms.RawTraces(:,iC))+2*std(ms.RawTraces(:,iC)));
    hold on
    plot(ms.time(t_binary)/1000, (ms.Binary(t_binary,iC)*0)+max(ms.RawTraces(:,iC)), '.', 'color', PARAMS.red)
    ylim([min(ylim), max(ylim)*1.2])
    
    
    %%% add the SPF for this cell.
    subplot(M, N, 3) % spf with centroid.
    Spr = winter(32);
    colormap([0 0 0 ; Spr(16:end,:)]);
    % c_lim = [0.2*max(max(ms.PeakToNoiseProj)), max(max(ms.PeakToNoiseProj))]; % helps clean up the projection by increasing the floor of the colormap to take in the lowest x% of the data
    % imagesc(ms.PeakToNoiseProj, c_lim)
    MS_plot_all_SFPs(flipdim(ms.SFPs,3)); % custom function to plot all the SFPs on top of each other.  Cleaner than ms.PeakToNoiseProj.
    hold on
    [max_I, max_J] = find(ms.SFPs(:,:,iC) == max(ms.SFPs(:,:,iC), [],[1,2]));
    text(max_J(1),max_I(1), '+', 'color', 'w',  'fontsize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    % quiver(max_J(1)-22,max_I(1)-3, 22,3,-10,'color', 'w', 'linewidth', 2, 'MaxHeadSize', 5); % add an arrow pointing to the current cell.
    %     scatter(ms.Centroids(iC,2), ms.Centroids(iC,1),60,'w', 'o','LineWidth',.5); % put a circle around the current cell.
    
    
    %%% place info
    % X Y position
    subplot(M, N, N+1:N+2)
    hold on
    plot(behav_aligned.time/1000, behav_aligned.position(:,1), 'color', PARAMS.L_grey)
    
    % update position in time with binary 'spikes'
    plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red);
    
    % plot(behav_aligned.time/1000, behav_aligned.position(:,2),'color', PARAMS.blue)
    ylabel('linear position')
    xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
    % legend({'x', 'y'})
    
    
    % plot the binary times on the position
    subplot(M, N, N+3) % N*4+4:N*4+6
    hold on
    plot(behav_aligned.position(:,1), behav_aligned.position(:,2), 'color', PARAMS.L_grey)
    set(gca, 'xtick', [], 'ytick', []);
    % put dots on positions when the cell was active.
    MS_color_plot(behav_aligned.position(t_binary,1), behav_aligned.position(t_binary,2), '.', cool(length(behav_aligned.position(t_binary,1))))
    axis off
 
    
    
    % add in the 2D place/spatial information?
    subplot(M, N, N+4) % N*4+4:N*4+6
    imagesc('XData', X_bin_centers, 'CData',Place_posterior');
    set(gca, 'YDir', 'normal');
    xlabel('position (cm)');
    ylabel('position (cm)');

    
    
    % plot the MI and p value for the cell.
    subplot(M, N, N+5)
    text(0, 1*max(ylim), 'Place', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
    text(0, .5*max(ylim), {'MI:'; num2str(All_cells.place.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    text(0, .1*max(ylim), {'split corr:'; num2str(All_cells.place.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    axis off
    
    
    
    %%% speed info
    subplot(M, N, N*2+1:N*2+2)
    hold on
    plot(behav_aligned.time/1000, behav_aligned.speed, 'color', PARAMS.L_grey)
    plot(behav_aligned.time(movement_idx)/1000, behav_aligned.speed(movement_idx),'.', 'color', PARAMS.gold, 'markersize', 1)
    % legend('Speed', 'box', 'off')
    
    xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
    ylabel('speed cm/s')
    xlabel('time (s)')
    
    % update speed in time with binary 'spikes'
    plot(behav_aligned.time(t_binary)/1000,behav_aligned.speed(t_binary,1),'.', 'color', PARAMS.red);
    
    % speed stats
    subplot(M, N, N*2+5)
    text(0, 1*max(ylim), 'Speed', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
    text(0, .5*max(ylim), {'MI:'; num2str(All_cells.speed.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    %     text(0, .4*max(ylim), {'split corr:'; num2str(All_cells.speed.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    axis off
    
    % p(active | speed)
    subplot(M, N, N*2+3:N*2+4)
    S_h = shadedErrorBar(Speed_bin_centers, All_cells.speed.stats{iC}.mean', [All_cells.speed.stats{iC}.upper_CI95'; All_cells.speed.stats{iC}.lower_CI95']);
    S_h.mainLine.Color = PARAMS.gold; S_h.mainLine.LineWidth = 2;
    
    S_h.patch.FaceAlpha = .6; % how transparent the shading will be.
    S_h.patch.FaceColor = PARAMS.L_grey;
    S_h.edge(1).Color = PARAMS.D_grey;
    S_h.edge(2).Color = PARAMS.D_grey;
    
    xlim([min(Speed_bin_centers) max(Speed_bin_centers)]);
    xlabel('speed (cm/s)');
    ylabel('P(act | speed)');
    
    %%% acceleration info
    subplot(M, N, N*3+1:N*3+2)
    hold on
    plot(behav_aligned.time(1:end-1)/1000, behav_aligned.accel, 'color', PARAMS.L_grey)
    plot(behav_aligned.time(accel_movement_idx)/1000, behav_aligned.accel(accel_movement_idx),'.', 'color', PARAMS.green, 'markersize', 1)
    
    xlim([behav_aligned.time(1)/1000 max(behav_aligned.time(1:end-1))/1000]);
    ylim([cfg.accel_bins(1) cfg.accel_bins(end)])
    ylabel('acceleration cm/s^2')
    xlabel('time (s)')
    
    %%% update accel in time with binary 'spikes'
    plot(behav_aligned.time(accel_t_binary)/1000,behav_aligned.accel(accel_t_binary,1),'.', 'color', PARAMS.red);
    
    % accel stats
    subplot(M, N, N*3+5)
    text(0, 1*max(ylim), 'Accel', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
    text(0, .5*max(ylim), {'MI:'; num2str(All_cells.accel.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    %     text(0, .2*max(ylim), {'split corr:'; num2str(All_cells.accel.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    axis off
    
    %plot the MI with CI
    subplot(M, N, N*3+3:N*3+4)
    A_h = shadedErrorBar(Acc_bin_centers, All_cells.accel.stats{iC}.mean', [All_cells.accel.stats{iC}.upper_CI95'; All_cells.accel.stats{iC}.lower_CI95']);
    A_h.mainLine.Color = PARAMS.green; A_h.mainLine.LineWidth = 2;
    
    A_h.patch.FaceAlpha = .6; % how transparent the shading will be.
    A_h.patch.FaceColor = PARAMS.L_grey;
    A_h.edge(1).Color = PARAMS.D_grey;
    A_h.edge(2).Color = PARAMS.D_grey;
    
    xlim([min(Acc_bin_centers) max(Acc_bin_centers)]);
    xlabel('acceleration (cm/s^2)');
    ylabel('P(act | acceleration)');
    
    % %%% orientation info
    % subplot(M, N, N*3+1:N*3+3)
    % plot(behav_aligned.time/1000,ones(size(behav_aligned.time)), 'color', 'w')
    % hold on
    % text(behav_aligned.time(floor(length(behav_aligned.time)/3))/1000, pi, 'HD placeholder')
    %     % ylabel('HD')
    % ylim([-pi pi])
    % set(gca, 'ytick', [-pi pi], 'yticklabel', {'-pi' 'pi'})
    
    
    
    % set(gca, 'yticklabel', num2str(roundn([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))],2)))
    % get the transient/position values
    % tran_x = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,1),ms.time(t_binary),'linear');
    % tran_y = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,2),ms.time(t_binary),'linear');
    %
    % plot(tran_x,tran_y,'.', 'color', PARAMS.red);
    
    % customize figure stuff
    
    pos = get(gcf, 'position');
    set(gcf, 'position', [pos(1)-pos(1)*.8 pos(2)-pos(2)*.8 pos(3)*2.7 pos(4) *1.8])
    tightfig
    
    % pause(3)
    %     close(100)
    [full,this_dir]=fileparts(pwd);
    [~,this_parent] = fileparts(full);
    mkdir([PARAMS.inter_dir filesep 'Summary' filesep this_parent filesep this_dir]);
%     saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep 'Cell_' num2str(iC) '_Spatial_info.fig'])
%     saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep 'Cell_' num2str(iC) '_Spatial_info.png'])
%     %
    
    %     close(300)
    fprintf('\n');
    
end % end cell loop.

%% make a plot of population level activity


% plot the population place scores.

figure(400)
hold on

%     bar(1:length(All_cells.place.MI))
subplot(3,1,1)
bar(1:length(All_cells.place.MI),All_cells.place.MI, 'facecolor', PARAMS.red)
%     xlabel(1:length(All_cells.place.MI));
ylabel('MI');

subplot(3,1,2)
bar(1:length(All_cells.place.MI),All_cells.place.split.Stability_corr,'facecolor', PARAMS.blue)
%     xlabel(1:length(All_cells.place.MI));
ylabel('Split corr');


subplot(3,1,3)
bar(1:length(All_cells.place.MI),All_cells.place.n_evts,'facecolor', PARAMS.gold)
xlabel('cell ID');
ylabel('n Transients');

    saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep  '_Pop_place.fig'])
    saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep  '_Pop_place.png'])

%population speed modulation

% bin pop binary signal into .5s bins (Chen 2015) to get % of cells active
% at a given time. get Peasron cor coeff between speed and pop activity.
% they normalized by mean and SD for both the speed vec and pop activity
% vec.
