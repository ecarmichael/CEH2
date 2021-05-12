function SI = Spatial_screener_info(cfg_in, f_info)
%% Spatial_screener_2D:
%
%
%
%    Inputs:
%    - cfg [struct]   configuration see the defaults below.
%
%    - f_info:  [struct]  contains experiment information, 
%
%
%
%    Outputs:
%    - SI [struct]  contains spatial information and stats for
%    place, speed, and acceleration.
%
%
%
%
% EC 2021-01-02   initial version
%
%% initialize

global PARAMS


cfg_def = [];
% general
cfg_def.method = 'Binary'; % can also be 'decon' for deconvolved. 
cfg_def.decon_thresh = .8; % threshold for deconvolved events. 
cfg_def.binary_thresh =2; % number of sd for binary threshold of zscored Ca transients.
cfg_def.CI = 0.95; % confidence intervals will be equal for lower as 1-cfg.CI
cfg_def.remove_cells = []; % cells to skip for analyses. NOT USED ATM
cfg_def.split_method = 'time'; % use time or nTrans for session splitting method.
cfg_def.units = 'cm'; %units for bin size. 
% place
cfg_def.p_thresh = 0.05; % value for pvalue cut off;
cfg_def.stability_thres = 0.5; % from van der Veldt 2020
cfg_def.nShuff = 1000;
cfg_def.p_bin_size = 3 ; % in cm


cfg_def.split_gaus_sd = 3; % sd for gaussian smoothing of place tuning for split session xcorr.

% speed
cfg_def.s_bin_size = 2.5;
cfg_def.s_bins  =  2.5:cfg_def.s_bin_size:25; %
cfg_def.s_bins(cfg_def.s_bins==0) = []; %remove 0 bin.

% acceleration
cfg_def.a_bin_size = .4;
cfg_def.a_bins  = -2:cfg_def.a_bin_size:2; % between -2cm/s^2 and 2cm/s^s with 20 bins matches van de Veldt et al. 2020
cfg_def.a_bins(cfg_def.a_bins==0) = []; %remove 0 bin.


cfg = ProcessConfig(cfg_def, cfg_in);

%% load data
if exist('behav_DLC.mat')
    load('behav_DLC.mat')
else
    load('behav.mat')
end
load('ms.mat')

% get the binary
ms = MS_msExtractBinary_detrendTraces(ms, cfg.binary_thresh);

% if using the deconvolved information replace the binary with it. 
if strcmpi(cfg.method, 'decon')
    addpath(PARAMS.OASIS_dir); 
    oasis_setup; % setup script
    
%     for iChan = length(ms.
    
end



figure(11)
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


%% get the bins
cfg.X_bins = 0:cfg.p_bin_size:ceil(max(behav_aligned.position(:,1)));
X_bin_centers = cfg.X_bins +  cfg.p_bin_size/2;
X_bin_centers = X_bin_centers(1:end-1);
% same for Y bins
cfg.Y_bins = 0:cfg.p_bin_size:ceil(max(behav_aligned.position(:,2)));
Y_bin_centers = cfg.Y_bins +  cfg.p_bin_size/2;
Y_bin_centers = Y_bin_centers(1:end-1);

% make speed bins
Speed_bin_centers = cfg.s_bins + cfg.s_bin_size/2;
Speed_bin_centers = Speed_bin_centers(1:end-1);

% acceleration bins
Acc_bin_centers = cfg.a_bins + cfg.a_bin_size/2;
Acc_bin_centers = Acc_bin_centers(1:end-1);

%% plot basics for each cell

for iC = 1:size(ms.Binary,2) % loop through cells
    fprintf('\nProcessing cell %d...    ', iC);
    

    % get basic spatial information
    
    % get place information
    
    if strcmpi(f_info.task, 'LT')
        [Place_MI, Place_posterior, Place_occupancy, Place_Prob_active, Place_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx,iC),behav_aligned.position(movement_idx,1), cfg.X_bins);
    else
        [Place_MI, Place_posterior, Place_occupancy, Place_Prob_active, Place_likelihood] = MS_get_spatial_information_2D(ms.Binary(movement_idx,iC),behav_aligned.position(movement_idx,:), cfg.X_bins, cfg.Y_bins );
    end
    
    % get speed information
    [Speed_MI, ~, Speed_occupancy,~, Speed_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx,iC), behav_aligned.speed(movement_idx), cfg.s_bins);
    
     % get acceleration information
    [Acc_MI, ~, Acc_occupancy,~, Acc_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1)), cfg.a_bins);

    
    %% compute shuffle values
    
    for iShuff   = cfg.nShuff:-1:1
        counter(iShuff, cfg.nShuff) % just a progress counter
        
        random_ts = ceil(rand*length(ms.time));
        shuffled_binarized = circshift(ms.Binary(movement_idx,iC),random_ts);
            if strcmpi(f_info.task, 'LT')
                [Place_Shuff_MI(iShuff), ~, ~, ~, Place_Shuff_likelihood(:,:,iShuff)] = MS_get_spatial_information(shuffled_binarized,behav_aligned.position(movement_idx,1), cfg.X_bins);
            else
                [Place_Shuff_MI(iShuff), ~, ~, ~, Place_Shuff_likelihood(:,:,iShuff)] = MS_get_spatial_information_2D(shuffled_binarized,behav_aligned.position(movement_idx,:), cfg.X_bins, cfg.Y_bins );
            end
        [Speed_Shuff_MI(iShuff), ~, ~,~, Speed_Shuff_likelihood(:,iShuff)] = MS_get_spatial_information(shuffled_binarized, behav_aligned.speed(movement_idx), cfg.s_bins);
        [Acc_Shuff_MI(iShuff), ~, ~,~, Acc_Shuff_likelihood(:,iShuff)] = MS_get_spatial_information(shuffled_binarized(1:end-1), behav_aligned.accel(movement_idx(1:end-1)), cfg.a_bins);
        
    end
    %% compute the significance values vs shuffle
     %get Place sig tuning
    Place_MI_pval = sum(Place_Shuff_MI > Place_MI,2)/cfg.nShuff;
    Place_MI_pval_z = (Place_MI -mean(Place_Shuff_MI))/std(Place_Shuff_MI);
    
    if strcmpi(f_info.task, 'LT')
        Place_map_pval = sum(Place_Shuff_likelihood > Place_likelihood,2)/cfg.nShuff;
    else
        Place_map_pval = sum(Place_Shuff_likelihood > Place_likelihood,3)/cfg.nShuff;
    end
    Place_Sig_map = Place_likelihood;
%     Place_Sig_map(Place_MI_pval > cfg.p_thresh) = 0;
    
    % Speed sig
    Speed_MI_pval = sum(Speed_Shuff_MI > Speed_MI,2)/cfg.nShuff;
    Speed_MI_pval_z = (Speed_MI -mean(Speed_Shuff_MI))/std(Speed_Shuff_MI);
    Speed_Sig_map = Speed_likelihood; 
%     Speed_Sig_map(Speed_MI_pval > cfg.p_thresh) = 0;
    
    
    % Acceleration sig
    Acc_MI_pval = sum(Acc_Shuff_MI > Acc_MI,2)/cfg.nShuff;
    Acc_MI_pval_z = (Acc_MI -mean(Acc_Shuff_MI))/std(Acc_Shuff_MI);
    Acc_Sig_map = Acc_likelihood;
%     Acc_Sig_map(Acc_MI_pval > cfg.p_thresh) = 0;
    %% compute the traditional (Solstad) border score WIP
    
%     MS_border_score(Place_Sig_map, 1, cfg.border_sd)

% %% compute EMD border score. 
    
       
    %% Compute the split half stats for place fields
%     fprintf('\n<strong>%s</strong>: computing split spatial info...\n', mfilename)
    % session splitting
    split_1 = zeros(size(ms.Binary(:,1)));
    ca_evts = MS_get_binary_events(ms.Binary(:,iC));
    
    if strcmp(cfg.split_method, 'ntran')
        % split based on number of transients;
        split_evt_idx = ca_evts(ceil(length(ca_evts)/2),1); % get the start of the middle event.
    else
        split_evt_idx = ceil(length(ms.Binary(:,1))/2); % hard split based on time.
    end
    
    split_1(1:split_evt_idx) = 1;
    
    % make keep indices for split halves.
    split_1 = logical(split_1);
    split_2 = logical(~split_1);
    
    % split 1
    if strcmpi(f_info.task, 'LT')
        [Place_S1_MI, ~, Place_S1_occupancy, ~, Place_S1_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx & split_1,iC),behav_aligned.position(movement_idx & split_1,1), cfg.X_bins);
        [Place_S1_shuff_MI, Place_S1_shuff_likelihood] = MS_split_shuff(ms.Binary(split_1,iC), behav_aligned.position(split_1,1),movement_idx, cfg.nShuff, cfg.X_bins);
    else
        [Place_S1_MI, ~, Place_S1_occupancy, ~, Place_S1_likelihood] = MS_get_spatial_information_2D(ms.Binary(movement_idx & split_1,iC),behav_aligned.position(movement_idx & split_1,:), cfg.X_bins, cfg.Y_bins );
        [Place_S1_shuff_MI, Place_S1_shuff_likelihood] = MS_split_shuff(ms.Binary(split_1,iC), behav_aligned.position(split_1,:),movement_idx, cfg.nShuff, cfg.X_bins, cfg.Y_bins);
    end
    
    
    Place_S1_MI_pval = sum(Place_S1_shuff_MI > Place_S1_MI,2)/cfg.nShuff;
    Place_S1_MI_pval_z = (Place_S1_MI -mean(Place_S1_shuff_MI))/std(Place_S1_shuff_MI);
    
    if strcmpi(f_info.task, 'LT')
        Place_S1_map_pval = sum(Place_S1_shuff_likelihood > Place_S1_likelihood,2)/cfg.nShuff;
    else
        Place_S1_map_pval = sum(Place_S1_shuff_likelihood > Place_S1_likelihood,3)/cfg.nShuff;
    end
    Place_S1_Sig_map = Place_S1_likelihood;
    Place_S1_Sig_map(Place_S1_MI_pval > cfg.p_thres) = 0;
    
    
    % split 2
   if strcmpi(f_info.task, 'LT')
        [Place_S2_MI, ~, Place_S2_occupancy, ~, Place_S2_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx & split_2,iC),behav_aligned.position(movement_idx & split_2,1), cfg.X_bins);
        [Place_S2_shuff_MI, Place_S2_shuff_likelihood] = MS_split_shuff(ms.Binary(split_2,iC), behav_aligned.position(split_2,1),movement_idx, cfg.nShuff, cfg.X_bins);
    else
        [Place_S2_MI, ~, Place_S2_occupancy, ~, Place_S2_likelihood] = MS_get_spatial_information_2D(ms.Binary(movement_idx & split_2,iC),behav_aligned.position(movement_idx & split_2,:), cfg.X_bins, cfg.Y_bins );
        [Place_S2_shuff_MI, Place_S2_shuff_likelihood] = MS_split_shuff(ms.Binary(split_2,iC), behav_aligned.position(split_2,:),movement_idx, cfg.nShuff, cfg.X_bins, cfg.Y_bins);
   end
    
   Place_S2_MI_pval = sum(Place_S2_shuff_MI > Place_S2_MI,2)/cfg.nShuff;
   Place_S2_MI_pval_z = (Place_S2_MI -mean(Place_S2_shuff_MI))/std(Place_S2_shuff_MI);
   
   if strcmpi(f_info.task, 'LT')
       Place_S2_map_pval = sum(Place_S2_shuff_likelihood > Place_S2_likelihood,2)/cfg.nShuff;
   else
       Place_S2_map_pval = sum(Place_S2_shuff_likelihood > Place_S2_likelihood,3)/cfg.nShuff;
   end
   
   Place_S2_Sig_map = Place_S2_likelihood;
   Place_S2_Sig_map(Place_S2_MI_pval > cfg.p_thres) = 0;
    
    
    % smooth with guassian
    S1_tuning_curve_smooth = imgaussfilt(Place_S1_likelihood, 2);
    S2_tuning_curve_smooth = imgaussfilt(Place_S2_likelihood, 2);
    Place_Stability_corr = corr2(S1_tuning_curve_smooth, S2_tuning_curve_smooth);
  
    
    %% collect information from each cell.
    SI(iC).finfo= f_info;
    SI(iC).cfg = cfg; 
    % ca events
    SI(iC).events.n_evts = length(ca_evts);
    SI(iC).events.ca_evts = ca_evts;
    
    % place info
    SI(iC).spatial.place.MI = Place_MI;
    SI(iC).spatial.place.MI_pval = Place_MI_pval;
    SI(iC).spatial.place.MI_pval_zscore = Place_MI_pval_z;
    SI(iC).spatial.place.Sig_map = Place_Sig_map;
    SI(iC).spatial.place.posterior = Place_posterior;
    SI(iC).spatial.place.probe_active = Place_Prob_active;
    SI(iC).spatial.place.occupany = Place_occupancy;
    SI(iC).spatial.place.likelihood = Place_likelihood;
    SI(iC).spatial.place.Shuff_MI = Place_Shuff_MI; 
    SI(iC).spatial.place.Shuff_lifelihood = Place_Shuff_likelihood; 

    
    %     SI.place.stats{iC} =  Place_stats;

    SI(iC).spatial.place.split.Stability_corr = Place_Stability_corr;
%     SI.place.split.split_sig(iC) = Place_split_sig;
%     SI.place.split.stats = Place_split_stats;
    
    SI(iC).spatial.place.split.S1_MI = Place_S1_MI;
    SI(iC).spatial.place.split.S1_Sig_map = Place_S1_Sig_map;
    SI(iC).spatial.place.split.S1_occupany = Place_S1_occupancy;
    
    SI(iC).spatial.place.split.S2_MI= Place_S2_MI;
    SI(iC).spatial.place.split.S2_Sig_map = Place_S2_Sig_map;
    SI(iC).spatial.place.split.S2_occupany = Place_S2_occupancy;
    
    
    % speed
    SI(iC).spatial.speed.MI = Speed_MI;
    SI(iC).spatial.speed.MI_pval = Speed_MI_pval;
    SI(iC).spatial.speed.MI_pval_zscore = Speed_MI_pval_z;
    SI(iC).spatial.speed.likelihood = Speed_likelihood;
    SI(iC).spatial.speed.Sig_map = Speed_Sig_map; 
    SI(iC).spatial.speed.occupany= Speed_occupancy;
%     SI.speed.stats{iC} = Speed_stats;
    
    % acceleration
    SI(iC).spatial.accel.MI = Acc_MI;
    SI(iC).spatial.accel.MI_pval = Acc_MI_pval;
    SI(iC).spatial.accel.MI_pval_zscore = Acc_MI_pval_z;
    SI(iC).spatial.accel.likelihood = Acc_likelihood;
    SI(iC).spatial.accel.Sig_map = Acc_Sig_map; 
    SI(iC).spatial.accel.occupany = Acc_occupancy;
%     SI.accel.stats{iC} = Acc_stats;
    
end
% save the output
save('SI.mat', 'SI');

%     
%     
%     %% figure 2 place information
%     %     if ishandle(200)
%     %         close(200)
%     %     end
%     figure(iC)
%     
%     subplot(6,4,[1 5])
%     axis off
%     text(0,1, ['Cell # ' num2str(iC)]);
%     text(0,.66,'Whole session')
%     text(0,.33,['MI: ' num2str(Place_MI,2)]);
%     text(0,0, ['Num trans: ' num2str(length(ca_evts))]);
%     colormap(gca, 'cool')
%     colorbar('location', 'southoutside', 'ticks', [0, 1], 'ticklabels', {'1^s^t', 'last'})
%     
%     
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,2)
%         axis off
%         title('Activity')
%         subplot(6,4,6)
%     else
%         subplot(6,4,[2 6])
%     end
%     t_binary = ms.Binary(:,iC) & movement_idx;
%     hold on
%     plot(behav_aligned.position(:,2), behav_aligned.position(:,1), 'color', PARAMS.L_grey)
%     %     xlim([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))])
%     %     ylim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]));
%     set(gca, 'xtick', [], 'ytick', []);
%     % put dots on positions when the cell was active.
%     MS_color_plot(behav_aligned.position(t_binary,2), behav_aligned.position(t_binary,1), '.', cool(length(behav_aligned.position(t_binary,2))))
%     %     plot(behav_aligned.position(t_binary,2), behav_aligned.position(t_binary,1),'.', 'color',  cool(length(behav_aligned.position(t_binary,2))))
%     %     xlim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]))
%     %     ylim(round([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))]))
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     hold on
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     % overall occupancy map
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,3)
%         title('Occupancy')
%         axis off
%         subplot(6,4,7);
%     else
%         subplot(6,4,[3 7]);
%         title('Occupancy')
%     end
%     imagesc(Y_bin_centers, X_bin_centers,  Place_occupancy);
%     set(gca,'YDir','normal'); % fix the Y direction to match the activity plot.
%     set(gca, 'xtick', [], 'ytick', []);
%     axis off
%     %     ylim([min(Y_bin_centers) max(Y_bin_centers)])
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     % overall tuning map
%     % overall occupancy map
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,4)
%         title(['Sig at p < ' num2str(cfg.p_thres)])
%         axis off
%         subplot(6,4,8);
%     else
%         subplot(6,4,[4 8]);
%         title(['Sig at p < ' num2str(cfg.p_thres)])
%     end
%     imagesc(Y_bin_centers,X_bin_centers, Place_Sig_map)
%     hold on
%     set(gca,'YDir','normal'); % fix the Y direction to match the activity plot.
%     set(gca, 'xtick', [], 'ytick', []);
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     
%     % first half
%     
%     subplot(6,4,[9 13])
%     axis off
%     text(0,.8,'1^s^t half split')
%     text(0,.6,['MI: ' num2str(Place_S1_MI,2)]);
%     
%     
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,14)
%     else
%         subplot(6,4,[10 14])
%     end
%     t_binary = ms.Binary(:,iC) & movement_idx & split_1;
%     hold on
%     plot(behav_aligned.position(movement_idx & split_1,2), behav_aligned.position(movement_idx & split_1,1), 'color', PARAMS.L_grey)
%     %     xlim([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))])
%     %     ylim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]));
%     set(gca, 'xtick', [], 'ytick', []);
%     % put dots on positions when the cell was active.
%     plot(behav_aligned.position(t_binary,2), behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red)
%     %     xlim(round([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))]))
%     %     ylim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]))
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     % overall occupancy map
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,15);
%     else
%         subplot(6,4,[11 15]);
%     end
%     imagesc(Y_bin_centers, X_bin_centers,  Place_S1_occupancy);
%     set(gca,'YDir','normal'); % fix the Y direction to match the activity plot
%     set(gca, 'xtick', [], 'ytick', []);
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     % 1st half tuning map
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,16);
%     else
%         subplot(6,4,[12 16]);
%     end
%     imagesc(Y_bin_centers, X_bin_centers,  Place_S2_likelihood);
%     set(gca,'YDir','normal'); % fix the Y direction to match the activity plot
%     set(gca, 'xtick', [], 'ytick', []);
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     
%     % second half
%     
%     subplot(6,4,17)
%     axis off
%     text(0,.8,'2^n^d half split')
%     text(0,.4,['MI: ' num2str(Place_S2_MI,2)]);
%     text(0,.0,['split xcorr: ' num2str(Place_Stability_corr,2)]);
%     
%     if contains(SI.fname{iC}.task, 'rec')
%         subplot(6,4,22)
%     else
%         subplot(6,4,[18 22])
%     end
%     t_binary = ms.Binary(:,iC) & movement_idx & split_2;
%     hold on
%     plot(behav_aligned.position(movement_idx & split_2,2), behav_aligned.position(movement_idx & split_2,1), 'color', PARAMS.L_grey)
%     %     xlim([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))])
%     %     ylim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]));
%     set(gca, 'xtick', [], 'ytick', []);
%     % put dots on positions when the cell was active.
%     plot(behav_aligned.position(t_binary,2), behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red)
%     %     xlim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]))
%     %     ylim(round([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))]))
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     % 2nd half occupancy
%     if contains(SI.fname{iC}.task, 'rec')
%         s2og = subplot(6,4,23);
%     else
%         s2og = subplot(6,4,[19 23]);
%     end
%     imagesc(Y_bin_centers,X_bin_centers, Place_S2_occupancy);
%     set(gca,'YDir','normal'); % fix the Y direction to match the activity plot
%     set(gca, 'xtick', [], 'ytick', []);
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     % overall tuning map
%     if contains(SI.fname{iC}.task, 'rec')
%         s2 =  subplot(6,4,24);
%     else
%         s2 =  subplot(6,4,[20 24]);
%     end
%     hold on
%     imagesc(X_bin_centers, Y_bin_centers, Place_S2_likelihood);
%     set(gca,'YDir','normal'); % fix the Y direction to match the activity plot
%     set(gca, 'xtick', [], 'ytick', []);
%     axis off
%     %     x_lim = xlim;
%     %     y_lim = ylim;
%     
%     %     set(gca, 'Clipping', 'off')
%     %     plot([x_lim(2)-8, x_lim(2)+2],[y_lim(1)-2, y_lim(1)-2], 'k', 'linewidth', 1)
%     %     plot([x_lim(2)+2, x_lim(2)+2],[y_lim(1)-2, y_lim(1)+8], 'k', 'linewidth', 1)
%     %     text(x_lim(2)-8, y_lim(1)-6, '10cm', 'fontsize', 6)
%     
%     %     % move unit line labels off
%     %     s1Pos = get(s2og, 'position');
%     %         s2Pos = get(s2, 'position');
%     %     s2Pos = [s2Pos(1)*.988 s2Pos(2) s1Pos(3:4)*1.4];
%     %     set(s2, 'position', s2Pos);
%     %     set(gcf, 'position', [680   421   886   550])
%     mkdir([PARAMS.inter_dir 'Place_figs'])
%     saveas(gcf, [PARAMS.inter_dir 'Place_figs' filesep f_info.subject '_' f_info.date '_' f_info.task '_Cell_' num2str(iC)], 'png')
%     saveas(gcf, [PARAMS.inter_dir  'Place_figs' filesep f_info.subject '_' f_info.date '_' f_info.task '_Cell_' num2str(iC)], 'fig')
%     %
%     %% plot everything
%     %     if ishandle(300)
%     %         close(300)
%     %     end
%     figure(100+iC)
%     
%     M = 4; % rows
%     N = 5; % columns
%     %   fig{iC} = figure('Visible', 'off'); % hack to stop figures from taking
%     %   over the screen.  Good for batch processing in the background.
%     
%     %get binary 'event times' to be plotted as dots
%     t_binary = ms.Binary(:,iC) & movement_idx;
%     %     accel_t_binary = find(ms.Binary(1:end-1,iC)==1);
%     accel_t_binary = t_binary(1:end-1);
%     
%     
%     %%% title information
%     subplot(M, N, 4:5) % title information. Top right corner,
%     ylim([0 10])
%     text(0,9,['Cell id: ' num2str(iC)], 'fontweight', 'bold')
%     text(0,7,['Subject: ' f_info.subject]);
%     text(0,5,['Session: ' f_info.task]);
%     text(0,3,['Date: ' f_info.date]);
%     text(0,1,['Binary thresh: ' num2str(ms.Binary_threshold) 'sd'  '    Num transients: ' num2str(length(ca_evts))])
%     
%     axis off
%     
%     
%     %%% raw trace
%     subplot(M, N, 1:2)
%     plot(ms.time/1000, ms.RawTraces(:,iC), 'color', PARAMS.blue)
%     xlim([ms.time(1)/1000 ms.time(end)/1000]);
%     xlabel('time(s)');
%     ylabel('dF/F');
%     hline(mean(ms.RawTraces(:,iC))+2*std(ms.RawTraces(:,iC)));
%     hold on
%     plot(ms.time(t_binary)/1000, (ms.Binary(t_binary,iC)*0)+max(ms.RawTraces(:,iC)), '.', 'color', PARAMS.red)
%     ylim([min(ylim), max(ylim)*1.2])
%     
%     
%     %%% add the SPF for this cell.
%     subplot(M, N, 3) % spf with centroid.
%     Spr = winter(32);
%     colormap([0 0 0 ; Spr(16:end,:)]);
%     % c_lim = [0.2*max(max(ms.PeakToNoiseProj)), max(max(ms.PeakToNoiseProj))]; % helps clean up the projection by increasing the floor of the colormap to take in the lowest x% of the data
%     % imagesc(ms.PeakToNoiseProj, c_lim)
%     MS_plot_all_SFPs(flipdim(ms.SFPs,3)); % custom function to plot all the SFPs on top of each other.  Cleaner than ms.PeakToNoiseProj.
%     hold on
%     [max_I, max_J] = find(ms.SFPs(:,:,iC) == max(ms.SFPs(:,:,iC), [],[1,2]));
%     text(max_J(1),max_I(1), '+', 'color', 'w',  'fontsize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
%     % quiver(max_J(1)-22,max_I(1)-3, 22,3,-10,'color', 'w', 'linewidth', 2, 'MaxHeadSize', 5); % add an arrow pointing to the current cell.
%     %     scatter(ms.Centroids(iC,2), ms.Centroids(iC,1),60,'w', 'o','LineWidth',.5); % put a circle around the current cell.
%     
%     
%     %%% place info
%     % X Y position
%     subplot(M, N, N+1:N+2)
%     hold on
%     plot(behav_aligned.time/1000, behav_aligned.position(:,1), 'color', PARAMS.L_grey)
%     plot(behav_aligned.time/1000, behav_aligned.position(:,2), 'color', PARAMS.D_grey)
%     
%     % update position in time with binary 'spikes'
%     plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red);
%     plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,2),'.', 'color', PARAMS.red);
%     
%     % plot(behav_aligned.time/1000, behav_aligned.position(:,2),'color', PARAMS.blue)
%     ylabel('linear position')
%     xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
%     % legend({'x', 'y'})
%     
%     
%     % plot the binary times on the position
%     subplot(M, N, N+3) % N*4+4:N*4+6
%     hold on
%     plot(behav_aligned.position(:,2), behav_aligned.position(:,1), 'color', PARAMS.L_grey)
%     set(gca, 'xtick', [], 'ytick', []);
%     % put dots on positions when the cell was active.
%     MS_color_plot(behav_aligned.position(t_binary,2), behav_aligned.position(t_binary,1), '.', cool(length(behav_aligned.position(t_binary,2))))
%     axis off
%     tmp=get(gca,'position'); % scale size of plot to match the
%     set(gca,'position',[tmp(1) tmp(2) tmp(3) (max(behav_aligned.position(:,1)) /max(behav_aligned.position(:,2)))*tmp(4)])
%     
%     
%     % add in the 2D place/spatial information?
%     subplot(M, N, N+4) % N*4+4:N*4+6
%     imagesc(Y_bin_centers,X_bin_centers,Place_posterior);
%     set(gca, 'YDir', 'normal');
%     xlabel('position (cm)');
%     ylabel('position (cm)');
%     tmp=get(gca,'position');
%     set(gca,'position',[tmp(1) tmp(2) tmp(3) (X_bin_centers(end) /Y_bin_centers(end))*tmp(4)])
%     
%     
%     % plot the MI and p value for the cell.
%     subplot(M, N, N+5)
%     text(0, 1*max(ylim), 'Place', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
%     text(0, .8*max(ylim), {'MI:'; num2str(SI.place.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
%     text(0, .4*max(ylim), {'split corr:'; num2str(SI.place.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
%     axis off
%     
%     
%     
%     %%% speed info
%     subplot(M, N, N*2+1:N*2+2)
%     hold on
%     plot(behav_aligned.time/1000, behav_aligned.speed, 'color', PARAMS.L_grey)
%     plot(behav_aligned.time(movement_idx)/1000, behav_aligned.speed(movement_idx),'.', 'color', PARAMS.gold, 'markersize', 1)
%     % legend('Speed', 'box', 'off')
%     
%     xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
%     ylabel('speed cm/s')
%     xlabel('time (s)')
%     
%     % update speed in time with binary 'spikes'
%     plot(behav_aligned.time(t_binary)/1000,behav_aligned.speed(t_binary,1),'.', 'color', PARAMS.red);
%     
%     % speed stats
%     subplot(M, N, N*2+5)
%     text(0, 1*max(ylim), 'Speed', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
%     text(0, .8*max(ylim), {'MI:'; num2str(SI.speed.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
%     %     text(0, .4*max(ylim), {'split corr:'; num2str(SI.speed.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
%     axis off
%     
%     % p(active | speed)
%     subplot(M, N, N*2+3:N*2+4)
%     S_h = shadedErrorBar(Speed_bin_centers, SI.speed.stats{iC}.mean', [SI.speed.stats{iC}.upper_CI95'; SI.speed.stats{iC}.lower_CI95']);
%     S_h.mainLine.Color = PARAMS.gold; S_h.mainLine.LineWidth = 2;
%     
%     S_h.patch.FaceAlpha = .6; % how transparent the shading will be.
%     S_h.patch.FaceColor = PARAMS.L_grey;
%     S_h.edge(1).Color = PARAMS.D_grey;
%     S_h.edge(2).Color = PARAMS.D_grey;
%     
%     xlim([min(Speed_bin_centers) max(Speed_bin_centers)]);
%     xlabel('speed (cm/s)');
%     ylabel('P(act | speed)');
%     
%     %%% acceleration info
%     subplot(M, N, N*3+1:N*3+2)
%     hold on
%     plot(behav_aligned.time(1:end-1)/1000, behav_aligned.accel, 'color', PARAMS.L_grey)
%     plot(behav_aligned.time(accel_movement_idx)/1000, behav_aligned.accel(accel_movement_idx),'.', 'color', PARAMS.green, 'markersize', 1)
%     
%     xlim([behav_aligned.time(1)/1000 max(behav_aligned.time(1:end-1))/1000]);
%     ylim([cfg.a_bins(1) cfg.a_bins(end)])
%     ylabel('acceleration cm/s^2')
%     xlabel('time (s)')
%     
%     %%% update accel in time with binary 'spikes'
%     plot(behav_aligned.time(accel_t_binary)/1000,behav_aligned.accel(accel_t_binary,1),'.', 'color', PARAMS.red);
%     
%     % accel stats
%     subplot(M, N, N*3+5)
%     text(0, 1*max(ylim), 'Accel', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
%     text(0, .6*max(ylim), {'MI:'; num2str(SI.accel.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
%     %     text(0, .2*max(ylim), {'split corr:'; num2str(SI.accel.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
%     axis off
%     
%     %plot the MI with CI
%     subplot(M, N, N*3+3:N*3+4)
%     A_h = shadedErrorBar(Acc_bin_centers, SI.accel.stats{iC}.mean', [SI.accel.stats{iC}.upper_CI95'; SI.accel.stats{iC}.lower_CI95']);
%     A_h.mainLine.Color = PARAMS.green; A_h.mainLine.LineWidth = 2;
%     
%     A_h.patch.FaceAlpha = .6; % how transparent the shading will be.
%     A_h.patch.FaceColor = PARAMS.L_grey;
%     A_h.edge(1).Color = PARAMS.D_grey;
%     A_h.edge(2).Color = PARAMS.D_grey;
%     
%     xlim([min(Acc_bin_centers) max(Acc_bin_centers)]);
%     xlabel('acceleration (cm/s^2)');
%     ylabel('P(act | acceleration)');
%     
%     % %%% orientation info
%     % subplot(M, N, N*3+1:N*3+3)
%     % plot(behav_aligned.time/1000,ones(size(behav_aligned.time)), 'color', 'w')
%     % hold on
%     % text(behav_aligned.time(floor(length(behav_aligned.time)/3))/1000, pi, 'HD placeholder')
%     %     % ylabel('HD')
%     % ylim([-pi pi])
%     % set(gca, 'ytick', [-pi pi], 'yticklabel', {'-pi' 'pi'})
%     
%     
%     
%     % set(gca, 'yticklabel', num2str(roundn([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))],2)))
%     % get the transient/position values
%     % tran_x = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,1),ms.time(t_binary),'linear');
%     % tran_y = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,2),ms.time(t_binary),'linear');
%     %
%     % plot(tran_x,tran_y,'.', 'color', PARAMS.red);
%     
%     % customize figure stuff
%     
%     pos = get(gcf, 'position');
%     set(gcf, 'position', [pos(1)-pos(1)*.8 pos(2)-pos(2)*.8 pos(3)*2.7 pos(4) *1.8])
%     tightfig
%     
%     % pause(3)
%     %     close(100)
%     [full,this_dir]=fileparts(pwd);
%     [~,this_parent] = fileparts(full);
%     
%     mkdir([PARAMS.inter_dir filesep 'Summary' filesep this_parent filesep this_dir]);
%     saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep 'Cell_' num2str(iC) '_Place_S1_shuff_likelihood.fig'])
%     saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep 'Cell_' num2str(iC) '_Place_S1_shuff_likelihood.png'])
%     %
%     
%     %     close(300)
%     fprintf('\n');
%     %     close(iC)
%     %     close(100+iC)
% end % end cell loop.
% 
% pause(1)
% %% make a plot of population level activity
% 
% 
% % plot the population place scores.
% 
% figure(400)
% hold on
% 
% %     bar(1:length(SI.place.MI))
% subplot(3,1,1)
% bar(1:length(SI.place.MI),SI.place.MI, 'facecolor', PARAMS.red)
% %     xlabel(1:length(SI.place.MI));
% ylabel('MI');
% 
% subplot(3,1,2)
% bar(1:length(SI.place.MI),SI.place.split.Stability_corr,'facecolor', PARAMS.blue)
% %     xlabel(1:length(SI.place.MI));
% ylabel('Split corr');
% 
% 
% subplot(3,1,3)
% bar(1:length(SI.place.MI),SI.place.n_evts,'facecolor', PARAMS.gold)
% xlabel('cell ID');
% ylabel('n Transients');
% 
% saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep  '_Pop_place.fig'])
% saveas(gcf, [PARAMS.inter_dir filesep 'Summary' filesep  this_parent filesep this_dir filesep  '_Pop_place.png'])
% 
% %population speed modulation
% 
% % bin pop binary signal into .5s bins (Chen 2015) to get % of cells active
% % at a given time. get Peasron cor coeff between speed and pop activity.
% % they normalized by mean and SD for both the speed vec and pop activity
% % vec.
