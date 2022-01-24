function [z_power, z_freq] = MS_zscore_LFP_JC(LFP_dir, CSC) 
%% MS_zscore_LFP_JC: loads the scored 'hypno' file, LFP, and miniscope LFP 
%  matrix, then computes the mean and SD from the filtered LFP signals 
%  before zscoring the values in the miniscope LFP matrices
%
%
%
%    Inputs: 
%    - none  [FUTURE put in specific valriables for flexibility] 
%
%
%
%    Outputs: 
%    - z_power:  [1 x nSamples array] 
%
%
%
%
% EC 2022-01-14   initial version 
%
%
%
%% initialize 
data_dir = cd; 
if exist(['pREM' filesep 'hypno.mat'], 'file')
    load(['pREM' filesep 'hypno.mat']);
else
    error('<strong>%s</strong>: no Hypno file detected in pREM folder\n', mfilename);
end
if ~exist(['LFP_mats' filesep 'all_d_post.mat'], 'file')
    error('<strong>%s</strong>: no LFP matrices detected in current directory\n', mfilename);
end





%% Load the CSC file and filter. 
parts = strsplit(cd,  filesep);

subject = parts{end-1};
type = strsplit(parts{end}, '_');
type = type{end};
if strcmp(type,'HATSwitch') && contains(subject, '1069')
    type = 'HATD6_switch';
elseif strcmp(type,'HATDSwitch') && contains(subject, '1060')
        type = 'HATSwitch';
end

% get the LFP files
cd(LFP_dir)

this_LFP_dir = MS_list_dir_names(cd, {subject, type});

cd(this_LFP_dir{1});

% if no CSC then load it using the same channel as the ms_seg. 
if isempty(CSC)
% load the MS_resize_seg to get tge the csc file
load([data_dir filesep 'ms_resize.mat']);

% get the channels to load from the pre-cut data. First channel should be
% EMG and second is the best LFP.
cfg_load = [];
for iC = length(ms_seg_resize.NLX_csc{1}.cfg.hdr) % loop over channels.
    cfg_load.fc{iC} = [ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.AcqEntName '.ncs'];
    cfg_load.desired_sampling_frequency = ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.SamplingFrequency;
end
cfg_load.fc(1) = []; 

% hard code LFP channel 
if strcmp(subject, 'PV1043') && strcmp(type, 'LTD5')
    cfg_load.fc{2} = 'CSC6.ncs';
end

% load some data.
CSC = MS_LoadCSC(cfg_load);

    Fs = floor(1/(0.001*(mode(diff((ms_seg_resize.time{1}))))));

else
    Fs = 30; % assume sampling rate for Ca data is 30fps. 
    
end
EVT = LoadEvents([]); % load the events file.

% clear ms_seg_resize

%% prepare the data for scoring.
% since JC data has a gap in a single recording for the track, account for
% that here.
S_rec_idx = find(contains(EVT.label, 'Starting Recording')); % get the index for the start recoding cell
Stop_rec_idx = find(contains(EVT.label, 'Stopping Recording')); % get the index for the start recoding cell

% split the pre and post recordings.  

if length(EVT.t{S_rec_idx}) ~= 2
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

 CSC_pre = CSC; 
    CSC_pre.tvec = CSC.tvec(1:nearest_idx(EVT.t{Stop_rec_idx}(pre_S_rec_idx),CSC.tvec));
    CSC_pre.data = CSC.data(:,1:nearest_idx(EVT.t{Stop_rec_idx}(pre_S_rec_idx),CSC.tvec)); 
%     CSC_pre.data = CSC_pre.data - mean(CSC_pre.data); 

    
    CSC_post= CSC;
    CSC_post.tvec = CSC.tvec(nearest_idx(EVT.t{S_rec_idx}(post_S_rec_idx), CSC.tvec):end);
    CSC_post.data = CSC.data(:,nearest_idx(EVT.t{S_rec_idx}(post_S_rec_idx), CSC.tvec):end);
%     CSC_post.data = CSC_post.data - mean(CSC_post.data);

CSC_cut = CSC;
CSC_cut.tvec = [CSC_pre.tvec; (CSC_post.tvec - CSC_post.tvec(1))+(CSC_pre.tvec(end)+(1/CSC.cfg.hdr{1}.SamplingFrequency))];
CSC_cut.data = [CSC_pre.data, CSC_post.data];

%% filter the CSC

  cfg_d = [];
    % filters
    cfg_d.type = 'fdesign'; %Cheby1 is sharper than butter
    cfg_d.f  = [1 5]; % broad
    cfg_d.order = 8; %type filter order (fine for this f range)
    cfg_d.display_filter = 0; % use this to see the fvtool
    delta = FilterLFP(cfg_d, CSC_cut);
    
    cfg_t = [];
    % filters
    cfg_t.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_t.f  = [6 11]; % 
    cfg_t.order = 3; %type filter order (fine for this f range)
    cfg_t.display_filter = 0; % use this to see the fvtool
    theta = FilterLFP(cfg_t, CSC_cut);
    
    cfg_lg = [];
    % filters
    cfg_lg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_lg.f  = [30 50]; % "mouse  low gamma'
    cfg_lg.order = 5; %type filter order (fine for this f range)
    cfg_lg.display_filter =0; % use this to see the fvtool
    LG = FilterLFP(cfg_lg, CSC_cut);
    
    cfg_mg = [];
    % filters
    cfg_mg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_mg.f  = [50 90]; % "mouse  low gamma'
    cfg_mg.order = 5; %type filter order (fine for this f range)
    cfg_mg.display_filter =0; % use this to see the fvtool
    MG = FilterLFP(cfg_mg, CSC_cut);
    
    cfg_hg = [];
    % filters
    cfg_hg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_hg.f  = [110 160]; %
    cfg_hg.order = 5; %type filter order (fine for this f range)
    cfg_hg.display_filter =0; % use this to see the fvtool
    HG = FilterLFP(cfg_hg, CSC_cut);
    
    cfg_uhg = [];
    % filters
    cfg_uhg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_uhg.f  = [160 250]; %
    cfg_uhg.order = 5; %type filter order (fine for this f range)
    cfg_uhg.display_filter =0; % use this to see the fvtool
    UHG = FilterLFP(cfg_uhg, CSC_cut);
    
    cfg_rip = [];
    % filters
    cfg_rip.type = 'butter'; %Cheby1 is sharper than butter
    cfg_rip.f  = [120 250]; % 
    cfg_rip.order = 4; %type filter order (fine for this f range)
    cfg_rip.display_filter =0; % use this to see the fvtool
    Ripple = FilterLFP(cfg_rip, CSC_cut);
    
    
    % smooth the amplitude using 100ms
    D_amp = smooth(abs(hilbert(delta.data)), floor(Fs*0.1));
    T_amp = smooth(abs(hilbert(theta.data)), floor(Fs*0.1));
    LG_amp = smooth(abs(hilbert(LG.data)), floor(Fs*0.1));
    MG_amp = smooth(abs(hilbert(MG.data)), floor(Fs*0.1));
    HG_amp = smooth(abs(hilbert(HG.data)), floor(Fs*0.1));
    UHG_amp = smooth(abs(hilbert(UHG.data)), floor(Fs*0.1));
    Rip_amp = smooth(abs(hilbert(Ripple.data)), floor(Fs*0.1));
    
    % get the frequency
    D_freq = MS_estimate_freq(delta.tvec, delta.data);
    T_freq = MS_estimate_freq(theta.tvec, theta.data);
    LG_freq = MS_estimate_freq(LG.tvec, LG.data);
    MG_freq = MS_estimate_freq(MG.tvec, MG.data);
    HG_freq = MS_estimate_freq(HG.tvec, HG.data);
    UHG_freq = MS_estimate_freq(UHG.tvec, UHG.data);
    Rip_freq = MS_estimate_freq(Ripple.tvec, Ripple.data);
    
    % zscore them 
    [~, D_amp_mu, D_amp_sig] = zscore(D_amp); 
    [~, T_amp_mu, T_amp_sig] = zscore(T_amp); 
    [~, LG_amp_mu, LG_amp_sig] = zscore(LG_amp); 
    [~, MG_amp_mu, MG_amp_sig] = zscore(MG_amp); 
    [~, HG_amp_mu, HG_amp_sig] = zscore(HG_amp); 
    [~, UHG_amp_mu, UHG_amp_sig] = zscore(UHG_amp); 
    [~, Rip_amp_mu, Rip_amp_sig] = zscore(Rip_amp); 
    
    % z frequency
    [~, D_freq_mu, D_freq_sig] = zscore(D_freq); 
    [~, T_freq_mu, T_freq_sig] = zscore(T_freq); 
    [~, LG_freq_mu, LG_freq_sig] = zscore(LG_freq); 
    [~, MG_freq_mu, MG_freq_sig] = zscore(MG_freq); 
    [~, HG_freq_mu, HG_freq_sig] = zscore(HG_freq); 
    [~, UHG_freq_mu, UHG_freq_sig] = zscore(UHG_freq); 
    [~, Rip_freq_mu, Rip_freq_sig] = zscore(Rip_freq); 
    
    
    % load, zscore, and copy the 
    cd(data_dir)
    mkdir('LFP_mats_z')
    
    %% all 
    %%%PRE
        % delta
    load(['LFP_mats' filesep 'all_d_pre.mat']);
    all_d_pre = (all_d_pre - D_amp_mu)./D_amp_sig;
    save(['LFP_mats_z' filesep 'all_d_pre.mat'], 'all_d_pre');
    
    load(['LFP_mats' filesep 'all_d_freq_pre.mat']);
    all_d_freq_pre = (all_d_freq_pre - D_freq_mu)./D_freq_sig;
    save(['LFP_mats_z' filesep 'all_d_freq_pre.mat'], 'all_d_freq_pre');
    
    % theta
    load(['LFP_mats' filesep 'all_t_pre.mat']);
    all_t_pre = (all_t_pre - T_amp_mu)./T_amp_sig;
    save(['LFP_mats_z' filesep 'all_t_pre.mat'], 'all_t_pre');
    
    load(['LFP_mats' filesep 'all_t_freq_pre.mat']);
    all_t_freq_pre = (all_t_freq_pre - T_freq_mu)./T_freq_sig;
    save(['LFP_mats_z' filesep 'all_t_freq_pre.mat'], 'all_t_freq_pre');
    
    % low gamma
    load(['LFP_mats' filesep 'all_LG_pre.mat']);
    all_LG_pre = (all_LG_pre - LG_amp_mu)./LG_amp_sig;
    save(['LFP_mats_z' filesep 'all_LG_pre.mat'], 'all_LG_pre');
    
    load(['LFP_mats' filesep 'all_LG_freq_pre.mat']);
    all_LG_freq_pre = (all_LG_freq_pre - LG_freq_mu)./LG_freq_sig;
    save(['LFP_mats_z' filesep 'all_LG_freq_pre.mat'], 'all_LG_freq_pre');
    
    % mid gamma
    load(['LFP_mats' filesep 'all_MG_pre.mat']);
    all_MG_pre = (all_MG_pre - MG_amp_mu)./MG_amp_sig;
    save(['LFP_mats_z' filesep 'all_MG_pre.mat'], 'all_MG_pre');
    
    load(['LFP_mats' filesep 'all_MG_freq_pre.mat']);
    all_MG_freq_pre = (all_MG_freq_pre - MG_freq_mu)./MG_freq_sig;
    save(['LFP_mats_z' filesep 'all_MG_freq_pre.mat'], 'all_MG_freq_pre');
    
    % high gamma
    load(['LFP_mats' filesep 'all_HG_pre.mat']);
    all_HG_pre = (all_HG_pre - HG_amp_mu)./HG_amp_sig;
    save(['LFP_mats_z' filesep 'all_HG_pre.mat'], 'all_HG_pre');
    
    load(['LFP_mats' filesep 'all_HG_freq_pre.mat']);
    all_HG_freq_pre = (all_HG_freq_pre - HG_freq_mu)./HG_freq_sig;
    save(['LFP_mats_z' filesep 'all_HG_freq_pre.mat'], 'all_HG_freq_pre');
    
    % ultra high gamma
    load(['LFP_mats' filesep 'all_UHG_pre.mat']);
    all_UHG_pre = (all_UHG_pre - UHG_amp_mu)./UHG_amp_sig;
    save(['LFP_mats_z' filesep 'all_UHG_pre.mat'], 'all_UHG_pre');
    
    load(['LFP_mats' filesep 'all_UHG_freq_pre.mat']);
    all_UHG_freq_pre = (all_UHG_freq_pre - UHG_freq_mu)./UHG_freq_sig;
    save(['LFP_mats_z' filesep 'all_UHG_freq_pre.mat'], 'all_UHG_freq_pre');
    
    % ripple
    load(['LFP_mats' filesep 'all_Rip_pre.mat']);
    all_Rip_pre = (all_Rip_pre - Rip_amp_mu)./Rip_amp_sig;
    save(['LFP_mats_z' filesep 'all_Rip_pre.mat'], 'all_Rip_pre');
    
    load(['LFP_mats' filesep 'all_Rip_freq_pre.mat']);
    all_Rip_freq_pre = (all_Rip_freq_pre - Rip_freq_mu)./Rip_freq_sig;
    save(['LFP_mats_z' filesep 'all_Rip_freq_pre.mat'], 'all_Rip_freq_pre');
    
    
    %%% POST
    % delta
    load(['LFP_mats' filesep 'all_d_post.mat']);
    all_d_post = (all_d_post - D_amp_mu)./D_amp_sig;
    save(['LFP_mats_z' filesep 'all_d_post.mat'], 'all_d_post');
    
    load(['LFP_mats' filesep 'all_d_freq_post.mat']);
    all_d_freq_post = (all_d_freq_post - D_freq_mu)./D_freq_sig;
    save(['LFP_mats_z' filesep 'all_d_freq_post.mat'], 'all_d_freq_post');
    
    % theta
    load(['LFP_mats' filesep 'all_t_post.mat']);
    all_t_post = (all_t_post - T_amp_mu)./T_amp_sig;
    save(['LFP_mats_z' filesep 'all_t_post.mat'], 'all_t_post');
    
    load(['LFP_mats' filesep 'all_t_freq_post.mat']);
    all_t_freq_post = (all_t_freq_post - T_freq_mu)./T_freq_sig;
    save(['LFP_mats_z' filesep 'all_t_freq_post.mat'], 'all_t_freq_post');
    
    % low gamma
    load(['LFP_mats' filesep 'all_LG_post.mat']);
    all_LG_post = (all_LG_post - LG_amp_mu)./LG_amp_sig;
    save(['LFP_mats_z' filesep 'all_LG_post.mat'], 'all_LG_post');
    
    load(['LFP_mats' filesep 'all_LG_freq_post.mat']);
    all_LG_freq_post = (all_LG_freq_post - LG_freq_mu)./LG_freq_sig;
    save(['LFP_mats_z' filesep 'all_LG_freq_post.mat'], 'all_LG_freq_post');
    
    % mid gamma
    load(['LFP_mats' filesep 'all_MG_post.mat']);
    all_MG_post = (all_MG_post - MG_amp_mu)./MG_amp_sig;
    save(['LFP_mats_z' filesep 'all_MG_post.mat'], 'all_MG_post');
    
    load(['LFP_mats' filesep 'all_MG_freq_post.mat']);
    all_MG_freq_post = (all_MG_freq_post - MG_freq_mu)./MG_freq_sig;
    save(['LFP_mats_z' filesep 'all_MG_freq_post.mat'], 'all_MG_freq_post');
    
    % high gamma
    load(['LFP_mats' filesep 'all_HG_post.mat']);
    all_HG_post = (all_HG_post - HG_amp_mu)./HG_amp_sig;
    save(['LFP_mats_z' filesep 'all_HG_post.mat'], 'all_HG_post');
    
    load(['LFP_mats' filesep 'all_HG_freq_post.mat']);
    all_HG_freq_post = (all_HG_freq_post - HG_freq_mu)./HG_freq_sig;
    save(['LFP_mats_z' filesep 'all_HG_freq_post.mat'], 'all_HG_freq_post');
    
    % ultra high gamma
    load(['LFP_mats' filesep 'all_UHG_post.mat']);
    all_UHG_post = (all_UHG_post - UHG_amp_mu)./UHG_amp_sig;
    save(['LFP_mats_z' filesep 'all_UHG_post.mat'], 'all_UHG_post');
    
    load(['LFP_mats' filesep 'all_UHG_freq_post.mat']);
    all_UHG_freq_post = (all_UHG_freq_post - UHG_freq_mu)./UHG_freq_sig;
    save(['LFP_mats_z' filesep 'all_UHG_freq_post.mat'], 'all_UHG_freq_post');
    
    % ripple
    load(['LFP_mats' filesep 'all_Rip_post.mat']);
    all_Rip_post = (all_Rip_post - Rip_amp_mu)./Rip_amp_sig;
    save(['LFP_mats_z' filesep 'all_Rip_post.mat'], 'all_Rip_post');
    
    load(['LFP_mats' filesep 'all_Rip_freq_post.mat']);
    all_Rip_freq_post = (all_Rip_freq_post - Rip_freq_mu)./Rip_freq_sig;
    save(['LFP_mats_z' filesep 'all_Rip_freq_post.mat'], 'all_Rip_freq_post');
    
    clear all*
    
    %% REM

    %%%PRE
        % delta
    load(['LFP_mats' filesep 'all_d_pre_REM.mat']);
    all_d_pre_REM = (all_d_pre_REM - D_amp_mu)./D_amp_sig;
    save(['LFP_mats_z' filesep 'all_d_pre_REM.mat'], 'all_d_pre_REM');
    
    load(['LFP_mats' filesep 'all_d_freq_pre_REM.mat']);
    all_d_freq_pre_REM = (all_d_freq_pre_REM - D_freq_mu)./D_freq_sig;
    save(['LFP_mats_z' filesep 'all_d_freq_pre_REM.mat'], 'all_d_freq_pre_REM');
    
    % theta
    load(['LFP_mats' filesep 'all_t_pre_REM.mat']);
    all_t_pre_REM = (all_t_pre_REM - T_amp_mu)./T_amp_sig;
    save(['LFP_mats_z' filesep 'all_t_pre_REM.mat'], 'all_t_pre_REM');
    
    load(['LFP_mats' filesep 'all_t_freq_pre_REM.mat']);
    all_t_freq_pre_REM = (all_t_freq_pre_REM - T_freq_mu)./T_freq_sig;
    save(['LFP_mats_z' filesep 'all_t_freq_pre_REM.mat'], 'all_t_freq_pre_REM');
    
    % low gamma
    load(['LFP_mats' filesep 'all_LG_pre_REM.mat']);
    all_LG_pre_REM = (all_LG_pre_REM - LG_amp_mu)./LG_amp_sig;
    save(['LFP_mats_z' filesep 'all_LG_pre_REM.mat'], 'all_LG_pre_REM');
    
    load(['LFP_mats' filesep 'all_LG_freq_pre_REM.mat']);
    all_LG_freq_pre_REM = (all_LG_freq_pre_REM - LG_freq_mu)./LG_freq_sig;
    save(['LFP_mats_z' filesep 'all_LG_freq_pre_REM.mat'], 'all_LG_freq_pre_REM');
    
    % mid gamma
    load(['LFP_mats' filesep 'all_MG_pre_REM.mat']);
    all_MG_pre_REM = (all_MG_pre_REM - MG_amp_mu)./MG_amp_sig;
    save(['LFP_mats_z' filesep 'all_MG_pre_REM.mat'], 'all_MG_pre_REM');
    
    load(['LFP_mats' filesep 'all_MG_freq_pre_REM.mat']);
    all_MG_freq_pre_REM = (all_MG_freq_pre_REM - MG_freq_mu)./MG_freq_sig;
    save(['LFP_mats_z' filesep 'all_MG_freq_pre_REM.mat'], 'all_MG_freq_pre_REM');
    
    % high gamma
    load(['LFP_mats' filesep 'all_HG_pre_REM.mat']);
    all_HG_pre_REM = (all_HG_pre_REM - HG_amp_mu)./HG_amp_sig;
    save(['LFP_mats_z' filesep 'all_HG_pre_REM.mat'], 'all_HG_pre_REM');
    
    load(['LFP_mats' filesep 'all_HG_freq_pre_REM.mat']);
    all_HG_freq_pre_REM = (all_HG_freq_pre_REM - HG_freq_mu)./HG_freq_sig;
    save(['LFP_mats_z' filesep 'all_HG_freq_pre_REM.mat'], 'all_HG_freq_pre_REM');
    
    % ultra high gamma
    load(['LFP_mats' filesep 'all_UHG_pre_REM.mat']);
    all_UHG_pre_REM = (all_UHG_pre_REM - UHG_amp_mu)./UHG_amp_sig;
    save(['LFP_mats_z' filesep 'all_UHG_pre_REM.mat'], 'all_UHG_pre_REM');
    
    load(['LFP_mats' filesep 'all_UHG_freq_pre_REM.mat']);
    all_UHG_freq_pre_REM = (all_UHG_freq_pre_REM - UHG_freq_mu)./UHG_freq_sig;
    save(['LFP_mats_z' filesep 'all_UHG_freq_pre_REM.mat'], 'all_UHG_freq_pre_REM');
    
    % ripple
    load(['LFP_mats' filesep 'all_Rip_pre_REM.mat']);
    all_Rip_pre_REM = (all_Rip_pre_REM - Rip_amp_mu)./Rip_amp_sig;
    save(['LFP_mats_z' filesep 'all_Rip_pre_REM.mat'], 'all_Rip_pre_REM');
    
    load(['LFP_mats' filesep 'all_Rip_freq_pre_REM.mat']);
    all_Rip_freq_pre_REM = (all_Rip_freq_pre_REM - Rip_freq_mu)./Rip_freq_sig;
    save(['LFP_mats_z' filesep 'all_Rip_freq_pre_REM.mat'], 'all_Rip_freq_pre_REM');
    
    
    %%% POST
    % delta
    load(['LFP_mats' filesep 'all_d_post_REM.mat']);
    all_d_post_REM = (all_d_post_REM - D_amp_mu)./D_amp_sig;
    save(['LFP_mats_z' filesep 'all_d_post_REM.mat'], 'all_d_post_REM');
    
    load(['LFP_mats' filesep 'all_d_freq_post_REM.mat']);
    all_d_freq_post_REM = (all_d_freq_post_REM - D_freq_mu)./D_freq_sig;
    save(['LFP_mats_z' filesep 'all_d_freq_post_REM.mat'], 'all_d_freq_post_REM');
    
    % theta
    load(['LFP_mats' filesep 'all_t_post_REM.mat']);
    all_t_post_REM = (all_t_post_REM - T_amp_mu)./T_amp_sig;
    save(['LFP_mats_z' filesep 'all_t_post_REM.mat'], 'all_t_post_REM');
    
    load(['LFP_mats' filesep 'all_t_freq_post_REM.mat']);
    all_t_freq_post_REM = (all_t_freq_post_REM - T_freq_mu)./T_freq_sig;
    save(['LFP_mats_z' filesep 'all_t_freq_post_REM.mat'], 'all_t_freq_post_REM');
    
    % low gamma
    load(['LFP_mats' filesep 'all_LG_post_REM.mat']);
    all_LG_post_REM = (all_LG_post_REM - LG_amp_mu)./LG_amp_sig;
    save(['LFP_mats_z' filesep 'all_LG_post_REM.mat'], 'all_LG_post_REM');
    
    load(['LFP_mats' filesep 'all_LG_freq_post_REM.mat']);
    all_LG_freq_post_REM = (all_LG_freq_post_REM - LG_freq_mu)./LG_freq_sig;
    save(['LFP_mats_z' filesep 'all_LG_freq_post_REM.mat'], 'all_LG_freq_post_REM');
    
    % mid gamma
    load(['LFP_mats' filesep 'all_MG_post_REM.mat']);
    all_MG_post_REM = (all_MG_post_REM - MG_amp_mu)./MG_amp_sig;
    save(['LFP_mats_z' filesep 'all_MG_post_REM.mat'], 'all_MG_post_REM');
    
    load(['LFP_mats' filesep 'all_MG_freq_post_REM.mat']);
    all_MG_freq_post_REM = (all_MG_freq_post_REM - MG_freq_mu)./MG_freq_sig;
    save(['LFP_mats_z' filesep 'all_MG_freq_post_REM.mat'], 'all_MG_freq_post_REM');
    
    % high gamma
    load(['LFP_mats' filesep 'all_HG_post_REM.mat']);
    all_HG_post_REM = (all_HG_post_REM - HG_amp_mu)./HG_amp_sig;
    save(['LFP_mats_z' filesep 'all_HG_post_REM.mat'], 'all_HG_post_REM');
    
    load(['LFP_mats' filesep 'all_HG_freq_post_REM.mat']);
    all_HG_freq_post_REM = (all_HG_freq_post_REM - HG_freq_mu)./HG_freq_sig;
    save(['LFP_mats_z' filesep 'all_HG_freq_post_REM.mat'], 'all_HG_freq_post_REM');
    
    % ultra high gamma
    load(['LFP_mats' filesep 'all_UHG_post_REM.mat']);
    all_UHG_post_REM = (all_UHG_post_REM - UHG_amp_mu)./UHG_amp_sig;
    save(['LFP_mats_z' filesep 'all_UHG_post_REM.mat'], 'all_UHG_post_REM');
    
    load(['LFP_mats' filesep 'all_UHG_freq_post_REM.mat']);
    all_UHG_freq_post_REM = (all_UHG_freq_post_REM - UHG_freq_mu)./UHG_freq_sig;
    save(['LFP_mats_z' filesep 'all_UHG_freq_post_REM.mat'], 'all_UHG_freq_post_REM');
    
    % ripple
    load(['LFP_mats' filesep 'all_Rip_post_REM.mat']);
    all_Rip_post_REM = (all_Rip_post_REM - Rip_amp_mu)./Rip_amp_sig;
    save(['LFP_mats_z' filesep 'all_Rip_post_REM.mat'], 'all_Rip_post_REM');
    
    load(['LFP_mats' filesep 'all_Rip_freq_post_REM.mat']);
    all_Rip_freq_post_REM = (all_Rip_freq_post_REM - Rip_freq_mu)./Rip_freq_sig;
    save(['LFP_mats_z' filesep 'all_Rip_freq_post_REM.mat'], 'all_Rip_freq_post_REM');
    
    %% SWS
    %%%PRE
        % delta
    load(['LFP_mats' filesep 'all_d_pre_SW.mat']);
    all_d_pre_SW = (all_d_pre_SW - D_amp_mu)./D_amp_sig;
    save(['LFP_mats_z' filesep 'all_d_pre_SW.mat'], 'all_d_pre_SW');
    
    load(['LFP_mats' filesep 'all_d_freq_pre_SW.mat']);
    all_d_freq_pre_SW = (all_d_freq_pre_SW - D_freq_mu)./D_freq_sig;
    save(['LFP_mats_z' filesep 'all_d_freq_pre_SW.mat'], 'all_d_freq_pre_SW');
    
    % theta
    load(['LFP_mats' filesep 'all_t_pre_SW.mat']);
    all_t_pre_SW = (all_t_pre_SW - T_amp_mu)./T_amp_sig;
    save(['LFP_mats_z' filesep 'all_t_pre_SW.mat'], 'all_t_pre_SW');
    
    load(['LFP_mats' filesep 'all_t_freq_pre_SW.mat']);
    all_t_freq_pre_SW = (all_t_freq_pre_SW - T_freq_mu)./T_freq_sig;
    save(['LFP_mats_z' filesep 'all_t_freq_pre_SW.mat'], 'all_t_freq_pre_SW');
    
    % low gamma
    load(['LFP_mats' filesep 'all_LG_pre_SW.mat']);
    all_LG_pre_SW = (all_LG_pre_SW - LG_amp_mu)./LG_amp_sig;
    save(['LFP_mats_z' filesep 'all_LG_pre_SW.mat'], 'all_LG_pre_SW');
    
    load(['LFP_mats' filesep 'all_LG_freq_pre_SW.mat']);
    all_LG_freq_pre_SW = (all_LG_freq_pre_SW - LG_freq_mu)./LG_freq_sig;
    save(['LFP_mats_z' filesep 'all_LG_freq_pre_SW.mat'], 'all_LG_freq_pre_SW');
    
    % mid gamma
    load(['LFP_mats' filesep 'all_MG_pre_SW.mat']);
    all_MG_pre_SW = (all_MG_pre_SW - MG_amp_mu)./MG_amp_sig;
    save(['LFP_mats_z' filesep 'all_MG_pre_SW.mat'], 'all_MG_pre_SW');
    
    load(['LFP_mats' filesep 'all_MG_freq_pre_SW.mat']);
    all_MG_freq_pre_SW = (all_MG_freq_pre_SW - MG_freq_mu)./MG_freq_sig;
    save(['LFP_mats_z' filesep 'all_MG_freq_pre_SW.mat'], 'all_MG_freq_pre_SW');
    
    % high gamma
    load(['LFP_mats' filesep 'all_HG_pre_SW.mat']);
    all_HG_pre_SW = (all_HG_pre_SW - HG_amp_mu)./HG_amp_sig;
    save(['LFP_mats_z' filesep 'all_HG_pre_SW.mat'], 'all_HG_pre_SW');
    
    load(['LFP_mats' filesep 'all_HG_freq_pre_SW.mat']);
    all_HG_freq_pre_SW = (all_HG_freq_pre_SW - HG_freq_mu)./HG_freq_sig;
    save(['LFP_mats_z' filesep 'all_HG_freq_pre_SW.mat'], 'all_HG_freq_pre_SW');
    
    % ultra high gamma
    load(['LFP_mats' filesep 'all_UHG_pre_SW.mat']);
    all_UHG_pre_SW = (all_UHG_pre_SW - UHG_amp_mu)./UHG_amp_sig;
    save(['LFP_mats_z' filesep 'all_UHG_pre_SW.mat'], 'all_UHG_pre_SW');
    
    load(['LFP_mats' filesep 'all_UHG_freq_pre_SW.mat']);
    all_UHG_freq_pre_SW = (all_UHG_freq_pre_SW - UHG_freq_mu)./UHG_freq_sig;
    save(['LFP_mats_z' filesep 'all_UHG_freq_pre_SW.mat'], 'all_UHG_freq_pre_SW');
    
    % ripple
    load(['LFP_mats' filesep 'all_Rip_pre_SW.mat']);
    all_Rip_pre_SW = (all_Rip_pre_SW - Rip_amp_mu)./Rip_amp_sig;
    save(['LFP_mats_z' filesep 'all_Rip_pre_SW.mat'], 'all_Rip_pre_SW');
    
    load(['LFP_mats' filesep 'all_Rip_freq_pre_SW.mat']);
    all_Rip_freq_pre_SW = (all_Rip_freq_pre_SW - Rip_freq_mu)./Rip_freq_sig;
    save(['LFP_mats_z' filesep 'all_Rip_freq_pre_SW.mat'], 'all_Rip_freq_pre_SW');
    
    
    %%% POST
    % delta
    load(['LFP_mats' filesep 'all_d_post_SW.mat']);
    all_d_post_SW = (all_d_post_SW - D_amp_mu)./D_amp_sig;
    save(['LFP_mats_z' filesep 'all_d_post_SW.mat'], 'all_d_post_SW');
    
    load(['LFP_mats' filesep 'all_d_freq_post_SW.mat']);
    all_d_freq_post_SW = (all_d_freq_post_SW - D_freq_mu)./D_freq_sig;
    save(['LFP_mats_z' filesep 'all_d_freq_post_SW.mat'], 'all_d_freq_post_SW');
    
    % theta
    load(['LFP_mats' filesep 'all_t_post_SW.mat']);
    all_t_post_SW = (all_t_post_SW - T_amp_mu)./T_amp_sig;
    save(['LFP_mats_z' filesep 'all_t_post_SW.mat'], 'all_t_post_SW');
    
    load(['LFP_mats' filesep 'all_t_freq_post_SW.mat']);
    all_t_freq_post_SW = (all_t_freq_post_SW - T_freq_mu)./T_freq_sig;
    save(['LFP_mats_z' filesep 'all_t_freq_post_SW.mat'], 'all_t_freq_post_SW');
    
    % low gamma
    load(['LFP_mats' filesep 'all_LG_post_SW.mat']);
    all_LG_post_SW = (all_LG_post_SW - LG_amp_mu)./LG_amp_sig;
    save(['LFP_mats_z' filesep 'all_LG_post_SW.mat'], 'all_LG_post_SW');
    
    load(['LFP_mats' filesep 'all_LG_freq_post_SW.mat']);
    all_LG_freq_post_SW = (all_LG_freq_post_SW - LG_freq_mu)./LG_freq_sig;
    save(['LFP_mats_z' filesep 'all_LG_freq_post_SW.mat'], 'all_LG_freq_post_SW');
    
    % mid gamma
    load(['LFP_mats' filesep 'all_MG_post_SW.mat']);
    all_MG_post_SW = (all_MG_post_SW - MG_amp_mu)./MG_amp_sig;
    save(['LFP_mats_z' filesep 'all_MG_post_SW.mat'], 'all_MG_post_SW');
    
    load(['LFP_mats' filesep 'all_MG_freq_post_SW.mat']);
    all_MG_freq_post_SW = (all_MG_freq_post_SW - MG_freq_mu)./MG_freq_sig;
    save(['LFP_mats_z' filesep 'all_MG_freq_post_SW.mat'], 'all_MG_freq_post_SW');
    
    % high gamma
    load(['LFP_mats' filesep 'all_HG_post_SW.mat']);
    all_HG_post_SW = (all_HG_post_SW - HG_amp_mu)./HG_amp_sig;
    save(['LFP_mats_z' filesep 'all_HG_post_SW.mat'], 'all_HG_post_SW');
    
    load(['LFP_mats' filesep 'all_HG_freq_post_SW.mat']);
    all_HG_freq_post_SW = (all_HG_freq_post_SW - HG_freq_mu)./HG_freq_sig;
    save(['LFP_mats_z' filesep 'all_HG_freq_post_SW.mat'], 'all_HG_freq_post_SW');
    
    % ultra high gamma
    load(['LFP_mats' filesep 'all_UHG_post_SW.mat']);
    all_UHG_post_SW = (all_UHG_post_SW - UHG_amp_mu)./UHG_amp_sig;
    save(['LFP_mats_z' filesep 'all_UHG_post_SW.mat'], 'all_UHG_post_SW');
    
    load(['LFP_mats' filesep 'all_UHG_freq_post_SW.mat']);
    all_UHG_freq_post_SW = (all_UHG_freq_post_SW - UHG_freq_mu)./UHG_freq_sig;
    save(['LFP_mats_z' filesep 'all_UHG_freq_post_SW.mat'], 'all_UHG_freq_post_SW');
    
    % ripple
    load(['LFP_mats' filesep 'all_Rip_post_SW.mat']);
    all_Rip_post_SW = (all_Rip_post_SW - Rip_amp_mu)./Rip_amp_sig;
    save(['LFP_mats_z' filesep 'all_Rip_post_SW.mat'], 'all_Rip_post_SW');
    
    load(['LFP_mats' filesep 'all_Rip_freq_post_SW.mat']);
    all_Rip_freq_post_SW = (all_Rip_freq_post_SW - Rip_freq_mu)./Rip_freq_sig;
    save(['LFP_mats_z' filesep 'all_Rip_freq_post_SW.mat'], 'all_Rip_freq_post_SW');
    