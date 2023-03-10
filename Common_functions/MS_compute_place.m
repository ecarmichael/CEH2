function data_out = MS_compute_place(cfg_in, ms, behav, f_info)
%% MS_compute_place: preprocesses the spatial properties of a recording session.
% Makes use of the simplified data extracted from the ms.mat and behav.mat from processed miniscope and DLC
%
%
%
%    Inputs:
%    - cfg_in: [struct]  contains configuration parameters. If fields are
%    specified in the cfg_in, they will overwrite the cfg_def defaults.
%
%    - ms: [struct]   simplified version of the ms data struct from CNFMe.
%
%    - behav: [struct]  output from DLC
%
%    - f_info: [struct] contains file information.
%
%    Outputs:
%    - data_out: [struct] contains ms, behav, behav_aligned, f_info, SI.
%
%
%
%
% EC 2021-08-18   initial version
% EC 2022-02-11   updates for W maze

%% initialize
cfg_def = [];
% general
cfg_def.cells = 1:size(ms.RawTraces,2); 
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
%% prepare data.

figure(11)
cfg_plot = [];
cfg_plot.view =[0 75];
cfg_plot.plot_type = '2d';
cfg_plot.Ca_type = 'detrendRaw';
cfg_plot.colors = parula((round(size(ms.RawTraces,2)/20)*1.1));
cfg_plot.Ca_chan = 1:round(size(ms.RawTraces,2)/20); % uncomment to plot all channels.
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

% if ~isfield(behav_aligned, 'speed')% uses same method as msExtractBehavior by GE;
    dt = nanmedian(diff(behav_aligned.time/1000));
    dx = [0; diff(behav_aligned.position(:,1))]; % GE modified
    dy = [0; diff(behav_aligned.position(:,2))]; % GE modified
    
    behav_aligned.speed = sqrt((dx).^2+(dy).^2)/dt;
    behav_aligned.speed(1) = behav_aligned.speed(2); % fill in fist idx with second due to diff.
% end
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

for iC = fliplr(cfg.cells) % loop through cells
    fprintf('\nProcessing cell %d...    ', iC);
    tic
    
    if sum(ms.Binary(movement_idx,iC)) == 0
        continue
    end
    % get basic spatial information
    
    % get place information
    [Place_MI, Place_posterior, Place_occupancy, Place_Prob_active, Place_likelihood] = MS_get_spatial_information_2D(ms.Binary(movement_idx,iC),behav_aligned.position(movement_idx,:), cfg.X_bins, cfg.Y_bins );
    
    % get speed information
    [Speed_MI, ~, Speed_occupancy,~, Speed_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx,iC), behav_aligned.speed(movement_idx), cfg.s_bins);
    
    % get acceleration information
    [Acc_MI, ~, Acc_occupancy,~, Acc_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx(1:end-1),iC), behav_aligned.accel(movement_idx(1:end-1)), cfg.a_bins);
    
    
    %% compute shuffle values
    % preallocate for parfor.
    task = f_info.task;
    nShuf = cfg.nShuff;
    nSamples = length(ms.time);
    
    Place_Shuff_MI  = NaN(1,cfg.nShuff); Speed_Shuff_MI  = NaN(1,cfg.nShuff); Acc_Shuff_MI  = NaN(1,cfg.nShuff);
  
        Place_Shuff_likelihood  = NaN(length(cfg.X_bins)-1, length(cfg.Y_bins)-1,cfg.nShuff);
    
    % shuffle
    tic
    fprintf('<strong>%s</strong>: shuffling...', mfilename);
    
    binary_mat = ms.Binary(movement_idx, iC); % for speed since parallel seems to take the entire struct. 
    pos_mat = behav_aligned.position(movement_idx,:); 
    spd_mat = behav_aligned.speed(movement_idx); 
    acc_mat =  behav_aligned.accel(movement_idx(1:end-1));
    parfor iShuff   = 1:nShuf
        %         counter(iShuff, nShuf) % just a progress counter
%         tic
%         disp(iShuff)
        random_ts = ceil(rand*nSamples);
        shuffled_binarized = circshift(binary_mat,random_ts);
        [Place_Shuff_MI(iShuff), ~, ~, ~, Place_Shuff_likelihood(:,:,iShuff)] = MS_get_spatial_information_2D(shuffled_binarized,pos_mat, cfg.X_bins, cfg.Y_bins );
        [Speed_Shuff_MI(iShuff), ~, ~,~, Speed_Shuff_likelihood(:,iShuff)] = MS_get_spatial_information(shuffled_binarized, spd_mat, cfg.s_bins);
        [Acc_Shuff_MI(iShuff), ~, ~,~, Acc_Shuff_likelihood(:,iShuff)] = MS_get_spatial_information(shuffled_binarized, acc_mat, cfg.a_bins);
%         toc
    end
    fprintf('\n')
    %% compute the significance values vs shuffle
    %get Place sig tuning
    Place_MI_pval = sum(Place_Shuff_MI > Place_MI,2)/cfg.nShuff;
    Place_MI_pval_z = (Place_MI -mean(Place_Shuff_MI))/std(Place_Shuff_MI);
    
    if strcmpi(f_info.task, 'LT') || strcmpi(f_info.task, 'Linear Track')
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
    if strcmpi(f_info.task, 'LT') || strcmpi(f_info.task, 'Linear Track')
        [Place_S1_MI, ~, Place_S1_occupancy, ~, Place_S1_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx & split_1,iC),behav_aligned.position(movement_idx & split_1,1), cfg.X_bins);
        [Place_S1_shuff_MI, Place_S1_shuff_likelihood] = MS_split_shuff(ms.Binary(split_1,iC), behav_aligned.position(split_1,1),movement_idx, cfg.nShuff, cfg.X_bins);
    else
        [Place_S1_MI, ~, Place_S1_occupancy, ~, Place_S1_likelihood] = MS_get_spatial_information_2D(ms.Binary(movement_idx & split_1,iC),behav_aligned.position(movement_idx & split_1,:), cfg.X_bins, cfg.Y_bins );
        [Place_S1_shuff_MI, Place_S1_shuff_likelihood] = MS_split_shuff(ms.Binary(split_1,iC), behav_aligned.position(split_1,:),movement_idx, cfg.nShuff, cfg.X_bins, cfg.Y_bins);
    end
    
    
    Place_S1_MI_pval = sum(Place_S1_shuff_MI > Place_S1_MI,2)/cfg.nShuff;
    Place_S1_MI_pval_z = (Place_S1_MI -mean(Place_S1_shuff_MI))/std(Place_S1_shuff_MI);
    
    if strcmpi(f_info.task, 'LT') || strcmpi(f_info.task, 'Linear Track')
        Place_S1_map_pval = sum(Place_S1_shuff_likelihood > Place_S1_likelihood,2)/cfg.nShuff;
    else
        Place_S1_map_pval = sum(Place_S1_shuff_likelihood > Place_S1_likelihood,3)/cfg.nShuff;
    end
    Place_S1_Sig_map = Place_S1_likelihood;
    Place_S1_Sig_map(Place_S1_MI_pval > cfg.p_thres) = 0;
    
    
    % split 2
    if strcmpi(f_info.task, 'LT') || strcmpi(f_info.task, 'Linear Track')
        [Place_S2_MI, ~, Place_S2_occupancy, ~, Place_S2_likelihood] = MS_get_spatial_information(ms.Binary(movement_idx & split_2,iC),behav_aligned.position(movement_idx & split_2,1), cfg.X_bins);
        [Place_S2_shuff_MI, Place_S2_shuff_likelihood] = MS_split_shuff(ms.Binary(split_2,iC), behav_aligned.position(split_2,1),movement_idx, cfg.nShuff, cfg.X_bins);
    else
        [Place_S2_MI, ~, Place_S2_occupancy, ~, Place_S2_likelihood] = MS_get_spatial_information_2D(ms.Binary(movement_idx & split_2,iC),behav_aligned.position(movement_idx & split_2,:), cfg.X_bins, cfg.Y_bins );
        [Place_S2_shuff_MI, Place_S2_shuff_likelihood] = MS_split_shuff(ms.Binary(split_2,iC), behav_aligned.position(split_2,:),movement_idx, cfg.nShuff, cfg.X_bins, cfg.Y_bins);
    end
    
    Place_S2_MI_pval = sum(Place_S2_shuff_MI > Place_S2_MI,2)/cfg.nShuff;
    Place_S2_MI_pval_z = (Place_S2_MI -mean(Place_S2_shuff_MI))/std(Place_S2_shuff_MI);
    
    if strcmpi(f_info.task, 'LT') || strcmpi(f_info.task, 'Linear Track')
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
    SI(iC).finfo.dir = pwd; % get the data dir.
    SI(iC).cfg = cfg;
    % ca events
    SI(iC).events.n_evts = length(ca_evts);
    SI(iC).events.ca_evts = ca_evts;
    
    % place info
    SI(iC).spatial.place.MI = Place_MI;
    SI(iC).spatial.place.MI_pval = Place_MI_pval;
    SI(iC).spatial.place.MI_pval_zscore = Place_MI_pval_z;
    SI(iC).spatial.place.Place_map_pval = Place_map_pval;

    SI(iC).spatial.place.Sig_map = Place_Sig_map';
    SI(iC).spatial.place.posterior = Place_posterior';
    SI(iC).spatial.place.probe_active = Place_Prob_active';
    SI(iC).spatial.place.occupany = Place_occupancy';
    SI(iC).spatial.place.likelihood = Place_likelihood';
    SI(iC).spatial.place.Shuff_MI = Place_Shuff_MI;
    SI(iC).spatial.place.Shuff_lifelihood = Place_Shuff_likelihood;
    
    %     SI.place.stats{iC} =  Place_stats;
    
    SI(iC).spatial.place.split.Stability_corr = Place_Stability_corr;
    %     SI.place.split.split_sig(iC) = Place_split_sig;
    %     SI.place.split.stats = Place_split_stats;
    
    SI(iC).spatial.place.split.S1_MI = Place_S1_MI;
    SI(iC).spatial.place.split.S1_Sig_map = Place_S1_Sig_map';
    SI(iC).spatial.place.split.S1_occupany = Place_S1_occupancy';
    SI(iC).spatial.place.split.S1_Place_S1_MI_pval_z = Place_S1_MI_pval_z;
    SI(iC).spatial.place.split.S1_Place_S1_map_pval = Place_S1_map_pval';

    SI(iC).spatial.place.split.S2_MI= Place_S2_MI;
    SI(iC).spatial.place.split.S2_Sig_map = Place_S2_Sig_map';
    SI(iC).spatial.place.split.S2_occupany = Place_S2_occupancy';
    SI(iC).spatial.place.split.S2_Place_S1_MI_pval_z = Place_S2_MI_pval_z;
    SI(iC).spatial.place.split.S2_Place_S1_map_pval = Place_S2_map_pval';

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
    toc
end

%% save the output.
data_out.f_info = f_info;
data_out.ms = ms;
data_out.behav = behav_aligned;
data_out.behav_aligned = behav_aligned;
data_out.SI = SI;
