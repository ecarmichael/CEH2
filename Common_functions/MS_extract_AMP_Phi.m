function ms_seg_resize = MS_extract_AMP_Phi(ms_load_dir, ms_save_dir, ms_fname_load, ms_fname_save,method)
%% MS_extract_AMP_Phi: repurposed MS_re_binarize_JC to only extract and conbine sleep LFP blocks.


%    Inputs:
%     - ms_load_dir [char]   path to ms_resize.mat file.
%
%     - ms_save_dir [char]   path to where you want the ms_resize.mat and the all_*.mat files to be updated and
%     saved.
%
%     - ms_fname_load [char] name of the .mat file to load.
%
%     - ms_fname_save [char] name of the .mat file to save
%
%     - method: [string] can be 'estimate' which will estimate the
%     frequency using the inter-peak interval.  Alternative is 'instfreq'
%     [WIP]
%
%    Outputs:
%     - ms_seg_resize: [struct] ms data that has been segmented.  contains
%     the updated ms_seg.resize.hypnolabel field.
%
%
%
%
% EC 2020-04-22   initial version
%
%
%% initialize

if nargin ==0
    ms_load_dir = cd;
    ms_save_dir = cd;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    method = 'estimate';
    fprintf('<strong>%s</strong>: no dir inputs, using cd: %s \n', mfilename, cd);
elseif nargin == 1
    ms_save_dir = ms_load_dir;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    method = 'estimate';
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin == 2
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    method = 'estimate';
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin ==3
    ms_fname_save = 'ms_resize'; % name of the file to save.
    method = 'estimate';
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
end

%% load data
if exist([ms_fname_load '.mat'], 'file')
    load(ms_fname_load);
else
    error('no ms_resize.mat found')
end

%% compute the power in different frequency bands if a csc is given as an input

for iSeg = length(ms_seg_resize.RawTraces):-1:1
    
    csc = ms_seg_resize.NLX_csc{iSeg};
    this_csc = csc;
    this_csc.data = csc.data(2,:);
    this_csc.label = [];
    this_csc.label{1} = csc.label{2};
    this_csc.cfg.hdr = [];
    this_csc.cfg.hdr{1} = csc.cfg.hdr{2};
    
    Fs = this_csc.cfg.hdr{1}.SamplingFrequency; 
    
    cfg_d = [];
    % filters
    cfg_d.type = 'fdesign'; %Cheby1 is sharper than butter
    cfg_d.f  = [1 5]; % broad
    cfg_d.order = 8; %type filter order (fine for this f range)
    cfg_d.display_filter = 0; % use this to see the fvtool
    delta = FilterLFP(cfg_d, this_csc);
    
    cfg_t = [];
    % filters
    cfg_t.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_t.f  = [5 12]; %
    cfg_t.order = 3; %type filter order (fine for this f range)
    cfg_t.display_filter = 0; % use this to see the fvtool
    theta = FilterLFP(cfg_t, this_csc);
    
    cfg_lg = [];
    % filters
    cfg_lg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_lg.f  = [30 55]; % "mouse  low gamma'
    cfg_lg.order = 5; %type filter order (fine for this f range)
    cfg_lg.display_filter =0; % use this to see the fvtool
    LG = FilterLFP(cfg_lg, this_csc);
    
    cfg_hg = [];
    % filters
    cfg_hg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_hg.f  = [70 90]; %
    cfg_hg.order = 5; %type filter order (fine for this f range)
    cfg_hg.display_filter =0; % use this to see the fvtool
    HG = FilterLFP(cfg_hg, this_csc);
    
    cfg_rip = [];
    % filters
    cfg_rip.type = 'butter'; %Cheby1 is sharper than butter
    cfg_rip.f  = [120 250]; %
    cfg_rip.order = 4; %type filter order (fine for this f range)
    cfg_rip.display_filter =0; % use this to see the fvtool
    Ripple = FilterLFP(cfg_rip, this_csc);
    
    Fs = floor(1/(0.001*(mode(diff((ms_seg_resize.time{1}))))));
    
    % smooth the amplitude using 100ms
    D_amp{iSeg} = smooth(abs(hilbert(delta.data)), floor(Fs*0.1));
    T_amp{iSeg} = smooth(abs(hilbert(theta.data)), floor(Fs*0.1));
    LG_amp{iSeg} = smooth(abs(hilbert(LG.data)), floor(Fs*0.1));
    HG_amp{iSeg} = smooth(abs(hilbert(HG.data)), floor(Fs*0.1));
    Rip_amp{iSeg} = smooth(abs(hilbert(Ripple.data)), floor(Fs*0.1));
    
    % get the frequency
%     if strcmp(method, 'estimate')
        D_freq{iSeg} = MS_estimate_freq(delta.tvec, delta.data);
        T_freq{iSeg} = MS_estimate_freq(theta.tvec, theta.data);
        LG_freq{iSeg} = MS_estimate_freq(LG.tvec, LG.data);
        HG_freq{iSeg} = MS_estimate_freq(HG.tvec, HG.data);
        Rip_freq{iSeg} = MS_estimate_freq(Ripple.tvec, Ripple.data);
%     else
%         D_freq{iSeg} = instfreq(delta.data,Fs,'FrequencyLimits',cfg_d.f(1):cfg_d.f(2));
%         
%         [p,f,t] = pspectrum(theta.data,Fs,'spectrogram','FrequencyLimits',cfg_t.f);
%         T_freq{iSeg} = instfreq(p,f, t);
%         LG_freq{iSeg} = instfreq(LG.data,Fs,'FrequencyLimits',cfg_lg.f(1):cfg_lg.f(2));
%         HG_freq{iSeg} = instfreq(HG.data,Fs,'FrequencyLimits',cfg_hg.f(1):cfg_hg.f(2));
%         Rip_freq{iSeg} = instfreq(Ripple.data,Fs,'FrequencyLimits',cfg_rip.f(1):cfg_rip.f(2));
%     end
    
    
    % 
    ms_seg = []; % cleared so that we can use this var name for saving.
    
    keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
    keep_idx =keep_idx(find((keep_idx ~= iSeg)));
    
    cfg_rem = [];
    ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
    
    % decell
    ms_seg = MS_de_cell(ms_seg);

    
    if length(ms_seg.time) ~= length(ms_seg.NLX_evt.t{end})
        fprintf('Segment # %s, length of time %d and NLX_events  %d do not match. Using ms_seg.time to generate NLX TS\n',  num2str(iSeg),length(ms_seg.time),length(ms_seg.NLX_evt.t{end}))
        c_ms_time = ((ms_seg.time-ms_seg.time(1))*0.001)';
        %         c_nlx_time = ((ms_seg.NLX_evt.t{end}-ms_seg.NLX_evt.t{end}(1)))%-0.0077;
        
        %         NLX_ts = [c_nlx_time, c_ms_time(length(c_nlx_time)+1:end)]+ms_seg.NLX_evt.t{end}(1);
        
        NLX_ts = c_ms_time + ms_seg.NLX_evt.t{end}(1);
    else
        NLX_ts = ms_seg.NLX_evt.t{end};
        
    end
    these_idx = nearest_idx3(NLX_ts,csc.tvec');
    
    % align the amplitude to the Ca time (frames)
    D_amp{iSeg} = D_amp{iSeg}(these_idx);
    T_amp{iSeg} = T_amp{iSeg}(these_idx);
    LG_amp{iSeg} = LG_amp{iSeg}(these_idx);
    HG_amp{iSeg} = HG_amp{iSeg}(these_idx);
    Rip_amp{iSeg} = Rip_amp{iSeg}(these_idx);
    
    % convert IPI to Freq.  
    D_freq{iSeg} = D_freq{iSeg}(these_idx);
    T_freq{iSeg} = T_freq{iSeg}(these_idx);
    LG_freq{iSeg} = LG_freq{iSeg}(these_idx);
    HG_freq{iSeg} = HG_freq{iSeg}(these_idx);
    Rip_freq{iSeg} = Rip_freq{iSeg}(these_idx);
    
end

%% recompute the all pre/post SW vs REM blocks.

all_d_pre = []; all_d_post = [];
all_t_pre = []; all_t_post = [];
all_LG_pre = []; all_LG_post = [];
all_HG_pre = []; all_HG_post = [];
all_Rip_pre = []; all_Rip_post = [];

all_d_pre_freq = []; all_d_post_freq = [];
all_t_pre_freq = []; all_t_post_freq = [];
all_LG_pre_freq = []; all_LG_post_freq = [];
all_HG_pre_freq = []; all_HG_post_freq = [];
all_Rip_pre_freq = []; all_Rip_post_freq = [];

% REM
all_d_pre_REM = []; all_d_post_REM = [];
all_t_pre_REM = []; all_t_post_REM = [];
all_LG_pre_REM = []; all_LG_post_REM = [];
all_HG_pre_REM = []; all_HG_post_REM = [];
all_Rip_pre_REM = []; all_Rip_post_REM = [];

all_d_pre_REM_freq = []; all_d_post_REM_freq = [];
all_t_pre_REM_freq = []; all_t_post_REM_freq = [];
all_LG_pre_REM_freq = []; all_LG_post_REM_freq = [];
all_HG_pre_REM_freq = []; all_HG_post_REM_freq = [];
all_Rip_pre_REM_freq = []; all_Rip_post_REM_freq = [];

% SWS
all_d_pre_SW = []; all_d_post_SW = [];
all_t_pre_SW = []; all_t_post_SW = [];
all_LG_pre_SW = []; all_LG_post_SW = [];
all_HG_pre_SW = []; all_HG_post_SW = [];
all_Rip_pre_SW = []; all_Rip_post_SW = [];

all_d_pre_SW_freq = []; all_d_post_SW_freq = [];
all_t_pre_SW_freq = []; all_t_post_SW_freq = [];
all_LG_pre_SW_freq = []; all_LG_post_SW_freq = [];
all_HG_pre_SW_freq = []; all_HG_post_SW_freq = [];
all_Rip_pre_SW_freq = []; all_Rip_post_SW_freq = [];


for iSeg = 1:length(ms_seg_resize.RawTraces)
    
    % cat the binary traces for pre V post, and REM v SW
    if strcmp(ms_seg_resize.pre_post{iSeg}, 'pre')
        % amplitude
        all_d_pre = [all_d_pre, D_amp{iSeg}'];
        all_t_pre = [all_t_pre, T_amp{iSeg}'];
        all_LG_pre = [all_LG_pre, LG_amp{iSeg}'];
        all_HG_pre = [all_HG_pre, HG_amp{iSeg}'];
        all_Rip_pre = [all_Rip_pre, Rip_amp{iSeg}'];
        
        % frequency
        all_d_pre_freq = [all_d_pre_freq, D_freq{iSeg}];
        all_t_pre_freq = [all_t_pre_freq, T_freq{iSeg}];
        all_LG_pre_freq = [all_LG_pre_freq, LG_freq{iSeg}];
        all_HG_pre_freq = [all_HG_pre_freq, HG_freq{iSeg}];
        all_Rip_pre_freq = [all_Rip_pre_freq, Rip_freq{iSeg}];
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            
            all_d_pre_REM = [all_d_pre_REM, D_amp{iSeg}'];
            all_t_pre_REM = [all_t_pre_REM, T_amp{iSeg}'];
            all_LG_pre_REM = [all_LG_pre_REM, LG_amp{iSeg}'];
            all_HG_pre_REM = [all_HG_pre_REM, HG_amp{iSeg}'];
            all_Rip_pre_REM = [all_Rip_pre_REM, Rip_amp{iSeg}'];
            
            % frequency
            all_d_pre_REM_freq = [all_d_pre_REM_freq, D_freq{iSeg}];
            all_t_pre_REM_freq = [all_t_pre_REM_freq, T_freq{iSeg}];
            all_LG_pre_REM_freq = [all_LG_pre_REM_freq, LG_freq{iSeg}];
            all_HG_pre_REM_freq = [all_HG_pre_REM_freq, HG_freq{iSeg}];
            all_Rip_pre_REM_freq = [all_Rip_pre_REM_freq, Rip_freq{iSeg}];
            
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            
            all_d_pre_SW = [all_d_pre_SW, D_amp{iSeg}'];
            all_t_pre_SW = [all_t_pre_SW, T_amp{iSeg}'];
            all_LG_pre_SW = [all_LG_pre_SW, LG_amp{iSeg}'];
            all_HG_pre_SW = [all_HG_pre_SW, HG_amp{iSeg}'];
            all_Rip_pre_SW = [all_Rip_pre_SW, Rip_amp{iSeg}'];
            
            % frequency
            all_d_pre_SW_freq = [all_d_pre_SW_freq, D_freq{iSeg}];
            all_t_pre_SW_freq = [all_t_pre_SW_freq, T_freq{iSeg}];
            all_LG_pre_SW_freq = [all_LG_pre_SW_freq, LG_freq{iSeg}];
            all_HG_pre_SW_freq = [all_HG_pre_SW_freq, HG_freq{iSeg}];
            all_Rip_pre_SW_freq = [all_Rip_pre_SW_freq, Rip_freq{iSeg}];
        end
        
    elseif strcmp(ms_seg_resize.pre_post{iSeg}, 'post')
        % amplitude
        all_d_post = [all_d_post, D_amp{iSeg}'];
        all_t_post = [all_t_post, T_amp{iSeg}'];
        all_LG_post = [all_LG_post, LG_amp{iSeg}'];
        all_HG_post = [all_HG_post, HG_amp{iSeg}'];
        all_Rip_post = [all_Rip_post, Rip_amp{iSeg}'];
        
        % frequency
        all_d_post_freq = [all_d_post_freq, D_freq{iSeg}];
        all_t_post_freq = [all_t_post_freq, T_freq{iSeg}];
        all_LG_post_freq = [all_LG_post_freq, LG_freq{iSeg}];
        all_HG_post_freq = [all_HG_post_freq, HG_freq{iSeg}];
        all_Rip_post_freq = [all_Rip_post_freq, Rip_freq{iSeg}];
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            
            all_d_post_REM = [all_d_post_REM, D_amp{iSeg}'];
            all_t_post_REM = [all_t_post_REM, T_amp{iSeg}'];
            all_LG_post_REM = [all_LG_post_REM, LG_amp{iSeg}'];
            all_HG_post_REM = [all_HG_post_REM, HG_amp{iSeg}'];
            all_Rip_post_REM = [all_Rip_post_REM, Rip_amp{iSeg}'];
            
            % frequency
            all_d_post_REM_freq = [all_d_post_REM_freq, D_freq{iSeg}];
            all_t_post_REM_freq = [all_t_post_REM_freq, T_freq{iSeg}];
            all_LG_post_REM_freq = [all_LG_post_REM_freq, LG_freq{iSeg}];
            all_HG_post_REM_freq = [all_HG_post_REM_freq, HG_freq{iSeg}];
            all_Rip_post_REM_freq = [all_Rip_post_REM_freq, Rip_freq{iSeg}];
            
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            
            all_d_post_SW = [all_d_post_SW, D_amp{iSeg}'];
            all_t_post_SW = [all_t_post_SW, T_amp{iSeg}'];
            all_LG_post_SW = [all_LG_post_SW, LG_amp{iSeg}'];
            all_HG_post_SW = [all_HG_post_SW, HG_amp{iSeg}'];
            all_Rip_post_SW = [all_Rip_post_SW, Rip_amp{iSeg}'];
            
            % frequency
            all_d_post_SW_freq = [all_d_post_SW_freq, D_freq{iSeg}];
            all_t_post_SW_freq = [all_t_post_SW_freq, T_freq{iSeg}];
            all_LG_post_SW_freq = [all_LG_post_SW_freq, LG_freq{iSeg}];
            all_HG_post_SW_freq = [all_HG_post_SW_freq, HG_freq{iSeg}];
            all_Rip_post_SW_freq = [all_Rip_post_SW_freq, Rip_freq{iSeg}];
        end
    end
end


%% save the files.
fprintf('<strong>%s</strong>: saving concatinated LFP amplitude and \n', mfilename);

% save the LFP arrays if needed

    lfp_mat_dir = [ms_save_dir filesep 'LFP_mats'];
    mkdir(lfp_mat_dir)
    f_list = {'d', 't', 'LG', 'HG', 'Rip'};
    
    for iF = 1:length(f_list)
        % pre
        %         save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre.mat'], ['all_' f_list{iF} '_pre'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_freq.mat'], ['all_' f_list{iF} '_pre_freq'], '-v7.3');
        
        
        % pre REM
        %         save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_REM.mat'], ['all_' f_list{iF} '_pre_REM'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_REM_freq.mat'], ['all_' f_list{iF} '_pre_REM_freq'], '-v7.3');
        
        % pre SW
        %         save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_SW.mat'], ['all_' f_list{iF} '_pre_SW'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_SW_freq.mat'], ['all_' f_list{iF} '_pre_SW_freq'], '-v7.3');
        
        
        %post
        %         save([lfp_mat_dir filesep 'all_' f_list{iF} '_post.mat'], ['all_' f_list{iF} '_post'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_freq.mat'], ['all_' f_list{iF} '_post_freq'], '-v7.3');
        
        % post REM
        %         save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_REM.mat'], ['all_' f_list{iF} '_post_REM'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_REM_freq.mat'], ['all_' f_list{iF} '_post_REM_freq'], '-v7.3');
        
        % post SW
        %         save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_SW.mat'], ['all_' f_list{iF} '_post_SW'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_SW_freq.mat'], ['all_' f_list{iF} '_post_SW_freq'], '-v7.3');
        
    end
    
%% save a file with the zscore config
% save([ms_save_dir filesep ms_fname_save '.mat'], 'ms_seg_resize', '-v7.3')
cfgs = [];
cfgs.date = datestr(date, 'yyyy_mm_dd');
cfgs.freq_method = method;

cfgs.filters.d = cfg_d;
cfgs.filters.t = cfg_t;
cfgs.filters.LG = cfg_lg;
cfgs.filters.HG = cfg_hg;
cfgs.filters.Rip = cfg_rip;

save([ms_save_dir filesep 'cfgs_lfp_on_' cfgs.date   '.mat'], 'cfgs', '-v7.3')

close all
end % end function