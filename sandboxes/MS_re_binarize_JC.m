%% re-score cut data sandbox

% %% the raw traces are here for each recording block:
% iRecord = 1;
% ms_seg_resize.RawTraces{iRecord}

% Binarized was only save back into the individual segments.  but we can
% make new ones here.

%Eric, before you make this script, i just wanto check,
%so this ms_seg is the session that we already cut and it has raw traces
%that is cut so, i guess when we change the binarize, we can just run for
%each of this, right?
% that is correct.  It would be very simple to just have a loop that loads,
% reruns the binary with a new threshold and then saves it back.  The one
% thing that would need to happen is to redo all the 'all_binary_..."
% parts.  which would be easy to do as well.
%Ok got it! i don't want to bother you for this, so i'll try to do it and
%upload my data! Or maybe guillaume might do it ? I just don't want to
%bother you haha

% writing the script would be very easy since it is just copying part of
% some other code and writing a quick loop.  One catch is that I can't do
% it tonight but I can have it done tomorrow.  it would not be any hassle.
%Ok! if you have time then i would be appreciate it but i can also do it if
%you busy anyway. i'll just upload the data today ! and that's it!  Okay
%sounds good.

%% get the file paths and loop

function ms_seg_resize = MS_re_binarize_JC(z_threshold, ms_load_dir, ms_save_dir, ms_fname_load, ms_fname_save, csc)
%% MS_reclass_hynpo: simple script for reclassifying the hypno labels in segemented ms data.
%
%
%
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
%     - csc can either be a real csc [struct] or just a string.  If struct
%     it will run the filter across the entire section and then extract the
%     samples that correspond to the miniscope frames.  Otherwise, it use
%     the NLX_csc field in the ms_resize.mat file which should contain the
%     raw LFP in the csc format.
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
    z_threshold = 3; % default from ms_binarizeDetrend
    ms_load_dir = cd;
    ms_save_dir = cd;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: no dir inputs, using cd: %s \n', mfilename, cd);
elseif nargin == 1
    ms_load_dir = cd;
    ms_save_dir = ms_load_dir;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin == 2
    ms_save_dir = ms_load_dir;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin ==3
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin == 4
    ms_fname_save = 'ms_resize.mat'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_save_dir);
end

%% load data
cd(ms_load_dir)
if exist([ms_fname_load '.mat'], 'file')
    load(ms_fname_load);
    
else
    error('no ms_resize.mat found')
end

%% load trk and update and extract time values.
if exist(['ms_trk.mat'], 'file')
    load('ms_trk.mat');
    
else
    error('no ms_trk.mat found')
end

% get time
time_str = regexp(ms_trk.file_names,'\d*','Match');
ms_trk.time_labels = datestr(datestr([time_str{1} ':' time_str{2} ':' time_str{3}]), 'HH:MM:SS');
temp_time = datetime(ms_trk.time_labels);
temp_time.Format = 'HH:mm:ss';
temp_time = temp_time+seconds((ms_trk.time(end)- ms_trk.time(1))/1000);

trk_time = datevec(ms_trk.time_labels);
trk_end_time = datevec(datestr(temp_time, 'HH:MM:SS'));

% save the trk segment back with the updated time field.
save('ms_trk.mat','ms_trk', '-v7.3')

% etime(trk_end_time, trk_time);
%% compute the power in different frequency bands if a csc is given as an input
if exist('csc','var') && isstruct(csc)
    
    cfg_d = [];
    % filters
    cfg_d.type = 'fdesign'; %Cheby1 is sharper than butter
    cfg_d.f  = [1 5]; % broad
    cfg_d.order = 8; %type filter order (fine for this f range)
    cfg_d.display_filter = 0; % use this to see the fvtool
    delta = FilterLFP(cfg_d, csc);
    
    cfg_t = [];
    % filters
    cfg_t.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_t.f  = [6 11]; % 
    cfg_t.order = 3; %type filter order (fine for this f range)
    cfg_t.display_filter = 0; % use this to see the fvtool
    theta = FilterLFP(cfg_t, csc);
    
    cfg_lg = [];
    % filters
    cfg_lg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_lg.f  = [30 50]; % "mouse  low gamma'
    cfg_lg.order = 5; %type filter order (fine for this f range)
    cfg_lg.display_filter =0; % use this to see the fvtool
    LG = FilterLFP(cfg_lg, csc);
    
    cfg_mg = [];
    % filters
    cfg_mg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_mg.f  = [50 90]; % "mouse  low gamma'
    cfg_mg.order = 5; %type filter order (fine for this f range)
    cfg_mg.display_filter =0; % use this to see the fvtool
    MG = FilterLFP(cfg_mg, csc);
    
    cfg_hg = [];
    % filters
    cfg_hg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_hg.f  = [110 160]; %
    cfg_hg.order = 5; %type filter order (fine for this f range)
    cfg_hg.display_filter =0; % use this to see the fvtool
    HG = FilterLFP(cfg_hg, csc);
    
    cfg_uhg = [];
    % filters
    cfg_uhg.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_uhg.f  = [160 250]; %
    cfg_uhg.order = 5; %type filter order (fine for this f range)
    cfg_uhg.display_filter =0; % use this to see the fvtool
    UHG = FilterLFP(cfg_uhg, csc);
    
    cfg_rip = [];
    % filters
    cfg_rip.type = 'butter'; %Cheby1 is sharper than butter
    cfg_rip.f  = [120 250]; % 
    cfg_rip.order = 4; %type filter order (fine for this f range)
    cfg_rip.display_filter =0; % use this to see the fvtool
    Ripple = FilterLFP(cfg_rip, csc);
    
    Fs = floor(1/(0.001*(mode(diff((ms_seg_resize.time{1}))))));
    
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
end

%% recompute the all pre/post SW vs REM blocks.
all_binary_pre = []; all_binary_post= [];
all_RawTraces_pre = []; all_RawTraces_post = [];
all_detrendRaw_pre = []; all_detrendRaw_post = [];
all_seg_idx = [];
%rem
all_binary_pre_REM = []; all_binary_post_REM= [];
all_RawTraces_pre_REM = []; all_RawTraces_post_REM = [];
all_detrendRaw_pre_REM = []; all_detrendRaw_post_REM = [];
pre_REM_idx = []; post_REM_idx = [];
% SW
all_binary_pre_SW = []; all_binary_post_SW= [];
all_RawTraces_pre_SW = []; all_RawTraces_post_SW = [];
all_detrendRaw_pre_SW = []; all_detrendRaw_post_SW = [];
pre_SW_idx = []; post_SW_idx = [];

if exist('csc', 'var')
    
    % AMPLITUDE
    all_d_pre = []; all_d_post = [];
    all_t_pre = []; all_t_post = [];
    all_LG_pre = []; all_LG_post = [];
    all_MG_pre = []; all_MG_post = [];
    all_HG_pre = []; all_HG_post = [];
    all_UHG_pre = []; all_UHG_post = [];
    all_Rip_pre = []; all_Rip_post = [];
    
    %     REM
    all_d_pre_REM = []; all_d_post_REM = [];
    all_t_pre_REM = []; all_t_post_REM = [];
    all_LG_pre_REM = []; all_LG_post_REM = [];
    all_MG_pre_REM = []; all_MG_post_REM = [];
    all_HG_pre_REM = []; all_HG_post_REM = [];
    all_UHG_pre_REM = []; all_UHG_post_REM = [];
    all_Rip_pre_REM = []; all_Rip_post_REM = [];
    
    % SWS
    all_d_pre_SW = []; all_d_post_SW = [];
    all_t_pre_SW = []; all_t_post_SW = [];
    all_LG_pre_SW = []; all_LG_post_SW = [];
    all_MG_pre_SW = []; all_MG_post_SW = [];
    all_HG_pre_SW = []; all_HG_post_SW = [];
    all_UHG_pre_SW = []; all_UHG_post_SW = [];
    all_Rip_pre_SW = []; all_Rip_post_SW = [];
    
    
    
    %%% FREQUENCY
    all_d_freq_pre = []; all_d_freq_post = [];
    all_t_freq_pre = []; all_t_freq_post = [];
    all_LG_freq_pre = []; all_LG_freq_post = [];
    all_MG_freq_pre = []; all_MG_freq_post = [];
    all_HG_freq_pre = []; all_HG_freq_post = [];
    all_UHG_freq_pre = []; all_UHG_freq_post = [];
    all_Rip_freq_pre = []; all_Rip_freq_post = [];
    
    %     REM
    all_d_freq_pre_REM = []; all_d_freq_post_REM = [];
    all_t_freq_pre_REM = []; all_t_freq_post_REM = [];
    all_LG_freq_pre_REM = []; all_LG_freq_post_REM = [];
    all_MG_freq_pre_REM = []; all_MG_freq_post_REM = [];
    all_HG_freq_pre_REM = []; all_HG_freq_post_REM = [];
    all_UHG_freq_pre_REM = []; all_UHG_freq_post_REM = [];
    all_Rip_freq_pre_REM = []; all_Rip_freq_post_REM = [];
    
    % SWS
    all_d_freq_pre_SW = []; all_d_freq_post_SW = [];
    all_t_freq_pre_SW = []; all_t_freq_post_SW = [];
    all_LG_freq_pre_SW = []; all_LG_freq_post_SW = [];
    all_MG_freq_pre_SW = []; all_MG_freq_post_SW = [];
    all_HG_freq_pre_SW = []; all_HG_freq_post_SW = [];
    all_UHG_freq_pre_SW = []; all_UHG_freq_post_SW = [];
    all_Rip_freq_pre_SW = []; all_Rip_freq_post_SW = [];
end


for iSeg = 1:length(ms_seg_resize.RawTraces)
    ms_seg = []; % cleared so that we can use this var name for saving.
    
    keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
    keep_idx =keep_idx(find((keep_idx ~= iSeg)));
    
    cfg_rem = [];
    ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
    
    ms_seg = MS_de_cell(ms_seg);
    
    % binarize the trace
    if ~isempty(z_threshold)
     ms_seg = MS_msExtractBinary_detrendTraces(ms_seg, z_threshold);
    end
    
   
    
    if length(ms_seg.time) ~= length(ms_seg.NLX_evt.t{end})
        fprintf('Segment # %s, length of time %d and NLX_events  %d do not match. Using ms_seg.time to generate NLX TS\n',  num2str(iSeg),length(ms_seg.time),length(ms_seg.NLX_evt.t{end}) )
        c_ms_time = ((ms_seg.time-ms_seg.time(1))*0.001)';
        %         c_nlx_time = ((ms_seg.NLX_evt.t{end}-ms_seg.NLX_evt.t{end}(1)))%-0.0077;
        
        %         NLX_ts = [c_nlx_time, c_ms_time(length(c_nlx_time)+1:end)]+ms_seg.NLX_evt.t{end}(1);
        
        NLX_ts = c_ms_time + ms_seg.NLX_evt.t{end}(1);
    else
        NLX_ts = ms_seg.NLX_evt.t{end};
        
    end
    these_idx = nearest_idx3(NLX_ts,csc.tvec');
    
    ms_seg.d_amp = D_amp(these_idx);
    ms_seg.t_amp = T_amp(these_idx);
    ms_seg.LG_amp = LG_amp(these_idx);
    ms_seg.MG_amp = MG_amp(these_idx);
    ms_seg.HG_amp = HG_amp(these_idx);
    ms_seg.UHG_amp = UHG_amp(these_idx);
    ms_seg.Rip_amp = Rip_amp(these_idx);
    
    % frequency
    ms_seg.d_freq = D_freq(these_idx);
    ms_seg.t_freq = T_freq(these_idx);
    ms_seg.LG_freq = LG_freq(these_idx);
    ms_seg.MG_freq = MG_freq(these_idx);
    ms_seg.HG_freq = HG_freq(these_idx);
    ms_seg.UHG_freq = UHG_freq(these_idx);
    ms_seg.Rip_freq = Rip_freq(these_idx);
    
    
    %     end
    
    %% to do add in thing that pulls out power for theta/lowG/highG/ripple?
    
    
    % add the spiking probability for different cell types.
    %     ms_seg.place_firing =
    
    % start of trk  =  trk_time; end of track = trk_end_time;
    % get the time vector for current segment.
    seg_time = datevec(ms_seg.time_labels);
    % add time vector
    if strcmp(ms_seg.pre_post, 'pre')
        ms_seg.time2trk = etime(seg_time, trk_time)/60;
    elseif strcmp(ms_seg.pre_post, 'post')
        ms_seg.time2trk = etime(seg_time,trk_end_time)/60;
    end
    
    cfg_SFP = [];
    cfg_SFP.fnc = '==';
    cfg_SFP.remove_val = 0;
    ms_seg = MS_update_SFP(cfg_SFP, ms_seg);
    
    this_dir = [];
    this_dir = [ms_save_dir filesep ms_seg_resize.file_names{iSeg}];
    fprintf('<strong>%s</strong>: saving resized ms struct back to %s...\n', mfilename, this_dir)
    if exist(this_dir, 'dir') % if the dir already exists, delete it and all the subfiles/dir and make a new one.  Avoids problems with multiple ms_resize.mat files in one dir.
        rmdir(this_dir, 's')
    end
    mkdir(this_dir)
    if  strcmp(ms_seg_resize.pre_post{iSeg}, 'post') && exist([this_dir filesep 'ms_seg_resize_' 'pre' '_' ms_seg_resize.hypno_label{iSeg} '.mat'])
        delete([this_dir filesep ms_fname_save '_' 'pre' '_' ms_seg_resize.hypno_label{iSeg} '.mat'])
    end
    % if this is a homecage do not use the 'pre' or post' lab.
    save([this_dir filesep ms_fname_save '_' ms_seg_resize.pre_post{iSeg} '_' ms_seg_resize.hypno_label{iSeg}],'ms_seg', '-v7.3');
    
    
    
    %% keep the index for the segment.
    if isempty(all_seg_idx)
        all_seg_idx(iSeg) = length(ms_seg.RawTraces);
    else
        all_seg_idx(iSeg) = length(ms_seg.RawTraces) + all_seg_idx(iSeg -1);
    end
    % cat the binary traces for pre V post, and REM v SW
    if strcmp(ms_seg_resize.pre_post{iSeg}, 'pre')
        all_binary_pre = [all_binary_pre; ms_seg.Binary];
        all_RawTraces_pre = [all_RawTraces_pre; ms_seg.RawTraces];
        all_detrendRaw_pre = [all_detrendRaw_pre; ms_seg.detrendRaw];
        
        if exist('csc', 'var')
            all_d_pre = [all_d_pre, ms_seg.d_amp'];
            all_t_pre = [all_t_pre, ms_seg.t_amp'];
            all_LG_pre = [all_LG_pre, ms_seg.LG_amp'];
            all_MG_pre = [all_MG_pre, ms_seg.MG_amp'];
            all_HG_pre = [all_HG_pre, ms_seg.HG_amp'];
            all_UHG_pre = [all_UHG_pre, ms_seg.UHG_amp'];
            all_Rip_pre = [all_Rip_pre, ms_seg.Rip_amp'];
            
            all_d_freq_pre = [all_d_freq_pre, ms_seg.d_freq];
            all_t_freq_pre = [all_t_freq_pre, ms_seg.t_freq];
            all_LG_freq_pre = [all_LG_freq_pre, ms_seg.LG_freq];
            all_MG_freq_pre = [all_MG_freq_pre, ms_seg.MG_freq];
            all_HG_freq_pre = [all_HG_freq_pre, ms_seg.HG_freq];
            all_UHG_freq_pre = [all_UHG_freq_pre, ms_seg.UHG_freq];
            all_Rip_freq_pre = [all_Rip_freq_pre, ms_seg.Rip_freq];
        end
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            all_binary_pre_REM = [all_binary_pre_REM; ms_seg.Binary];
            all_RawTraces_pre_REM = [all_RawTraces_pre_REM; ms_seg.RawTraces];
            all_detrendRaw_pre_REM = [all_detrendRaw_pre_REM; ms_seg.detrendRaw];
            pre_REM_idx = [pre_REM_idx, iSeg];
            
            if exist('csc', 'var')
                all_d_pre_REM = [all_d_pre_REM, ms_seg.d_amp'];
                all_t_pre_REM = [all_t_pre_REM, ms_seg.t_amp'];
                all_LG_pre_REM = [all_LG_pre_REM, ms_seg.LG_amp'];
                all_MG_pre_REM = [all_MG_pre_REM, ms_seg.MG_amp'];
                all_HG_pre_REM = [all_HG_pre_REM, ms_seg.HG_amp'];
                all_UHG_pre_REM = [all_UHG_pre_REM, ms_seg.UHG_amp'];
                all_Rip_pre_REM = [all_Rip_pre, ms_seg.Rip_amp'];
                
                all_d_freq_pre_REM = [all_d_freq_pre_REM, ms_seg.d_freq];
                all_t_freq_pre_REM = [all_t_freq_pre_REM, ms_seg.t_freq];
                all_LG_freq_pre_REM = [all_LG_freq_pre_REM, ms_seg.LG_freq];
                all_MG_freq_pre_REM = [all_MG_freq_pre_REM, ms_seg.MG_freq];
                all_HG_freq_pre_REM = [all_HG_freq_pre_REM, ms_seg.HG_freq];
                all_UHG_freq_pre_REM = [all_UHG_freq_pre_REM, ms_seg.UHG_freq];
                all_Rip_freq_pre_REM = [all_Rip_freq_pre, ms_seg.Rip_freq];
            end
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            all_binary_pre_SW = [all_binary_pre_SW; ms_seg.Binary];
            all_RawTraces_pre_SW = [all_RawTraces_pre_SW; ms_seg.RawTraces];
            all_detrendRaw_pre_SW = [all_detrendRaw_pre_SW; ms_seg.detrendRaw];
            pre_SW_idx = [pre_SW_idx, iSeg];
            if exist('csc', 'var')
                
                all_d_pre_SW = [all_d_pre_SW, ms_seg.d_amp'];
                all_t_pre_SW = [all_t_pre_SW, ms_seg.t_amp'];
                all_LG_pre_SW = [all_LG_pre_SW, ms_seg.LG_amp'];
                all_MG_pre_SW = [all_MG_pre_SW, ms_seg.MG_amp'];
                all_HG_pre_SW = [all_HG_pre_SW, ms_seg.HG_amp'];
                all_UHG_pre_SW = [all_UHG_pre_SW, ms_seg.UHG_amp'];
                all_Rip_pre_SW = [all_Rip_pre_SW, ms_seg.Rip_amp'];
                
                all_d_freq_pre_SW = [all_d_freq_pre_SW, ms_seg.d_freq];
                all_t_freq_pre_SW = [all_t_freq_pre_SW, ms_seg.t_freq];
                all_LG_freq_pre_SW = [all_LG_freq_pre_SW, ms_seg.LG_freq];
                all_MG_freq_pre_SW = [all_MG_freq_pre_SW, ms_seg.MG_freq];
                all_HG_freq_pre_SW = [all_HG_freq_pre_SW, ms_seg.HG_freq];
                all_UHG_freq_pre_SW = [all_UHG_freq_pre_SW, ms_seg.UHG_freq];
                all_Rip_freq_pre_SW = [all_Rip_freq_pre_SW, ms_seg.Rip_freq];
            end
        end
        
    elseif strcmp(ms_seg_resize.pre_post{iSeg}, 'post')
        all_binary_post = [all_binary_post; ms_seg.Binary];
        all_RawTraces_post = [all_RawTraces_post; ms_seg.RawTraces];
        all_detrendRaw_post = [all_detrendRaw_post; ms_seg.detrendRaw];
        if exist('csc', 'var')
            all_d_post = [all_d_post, ms_seg.d_amp'];
            all_t_post = [all_t_post, ms_seg.t_amp'];
            all_LG_post = [all_LG_post, ms_seg.LG_amp'];
            all_MG_post = [all_MG_post, ms_seg.MG_amp'];    
            all_HG_post = [all_HG_post, ms_seg.HG_amp'];
            all_UHG_post = [all_UHG_post, ms_seg.UHG_amp'];
            all_Rip_post = [all_Rip_post, ms_seg.Rip_amp'];
            
            all_d_freq_post = [all_d_freq_post, ms_seg.d_freq];
            all_t_freq_post = [all_t_freq_post, ms_seg.t_freq];
            all_LG_freq_post = [all_LG_freq_post, ms_seg.LG_freq];
            all_MG_freq_post = [all_MG_freq_post, ms_seg.MG_freq];    
            all_HG_freq_post = [all_HG_freq_post, ms_seg.HG_freq];
            all_UHG_freq_post = [all_UHG_freq_post, ms_seg.UHG_freq];
            all_Rip_freq_post = [all_Rip_freq_post, ms_seg.Rip_freq];
        end
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            all_binary_post_REM = [all_binary_post_REM; ms_seg.Binary];
            all_RawTraces_post_REM = [all_RawTraces_post_REM; ms_seg.RawTraces];
            all_detrendRaw_post_REM = [all_detrendRaw_post_REM; ms_seg.detrendRaw];
            post_REM_idx = [post_REM_idx, iSeg];
            if exist('csc', 'var')
                all_d_post_REM = [all_d_post_REM, ms_seg.d_amp'];
                all_t_post_REM = [all_t_post_REM, ms_seg.t_amp'];
                all_LG_post_REM = [all_LG_post_REM, ms_seg.LG_amp'];
                all_MG_post_REM = [all_MG_post_REM, ms_seg.MG_amp'];
                all_HG_post_REM = [all_HG_post_REM, ms_seg.HG_amp'];
                all_UHG_post_REM = [all_UHG_post_REM, ms_seg.UHG_amp'];
                all_Rip_post_REM = [all_Rip_post_REM, ms_seg.Rip_amp'];
                
                all_d_freq_post_REM = [all_d_freq_post_REM, ms_seg.d_freq];
                all_t_freq_post_REM = [all_t_freq_post_REM, ms_seg.t_freq];
                all_LG_freq_post_REM = [all_LG_freq_post_REM, ms_seg.LG_freq];
                all_MG_freq_post_REM = [all_MG_freq_post_REM, ms_seg.MG_freq];
                all_HG_freq_post_REM = [all_HG_freq_post_REM, ms_seg.HG_freq];
                all_UHG_freq_post_REM = [all_UHG_freq_post_REM, ms_seg.UHG_freq];
                all_Rip_freq_post_REM = [all_Rip_freq_post_REM, ms_seg.Rip_freq];
            end
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            all_binary_post_SW = [all_binary_post_SW; ms_seg.Binary];
            all_RawTraces_post_SW = [all_RawTraces_post_SW; ms_seg.RawTraces];
            all_detrendRaw_post_SW = [all_detrendRaw_post_SW; ms_seg.detrendRaw];
            post_SW_idx = [post_SW_idx, iSeg];
            if exist('csc', 'var')
                all_d_post_SW = [all_d_post_SW, ms_seg.d_amp'];
                all_t_post_SW = [all_t_post_SW, ms_seg.t_amp'];
                all_LG_post_SW = [all_LG_post_SW, ms_seg.LG_amp'];
                all_MG_post_SW = [all_MG_post_SW, ms_seg.MG_amp'];
                all_HG_post_SW = [all_HG_post_SW, ms_seg.HG_amp'];
                all_UHG_post_SW = [all_UHG_post_SW, ms_seg.UHG_amp'];
                all_Rip_post_SW = [all_Rip_post_SW, ms_seg.Rip_amp'];
                
                all_d_freq_post_SW = [all_d_freq_post_SW, ms_seg.d_freq];
                all_t_freq_post_SW = [all_t_freq_post_SW, ms_seg.t_freq];
                all_LG_freq_post_SW = [all_LG_freq_post_SW, ms_seg.LG_freq];
                all_MG_freq_post_SW = [all_MG_freq_post_SW, ms_seg.MG_freq];
                all_HG_freq_post_SW = [all_HG_freq_post_SW, ms_seg.HG_freq];
                all_UHG_freq_post_SW = [all_UHG_freq_post_SW, ms_seg.UHG_freq];
                all_Rip_freq_post_SW = [all_Rip_freq_post_SW, ms_seg.Rip_freq];
            end
        end
    end
end
all_seg_idx = [0 all_seg_idx];


%% save the files.
fprintf('<strong>%s</strong>: saving concatinating Binary, RawTraces, detrendRaw, and indicies\n', mfilename);


% save everything
save([ms_save_dir filesep 'all_seg_idx.mat'], 'all_seg_idx', '-v7.3');

% pre only
save([ms_save_dir filesep 'all_binary_pre.mat'], 'all_binary_pre', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_pre.mat'], 'all_RawTraces_pre', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_pre.mat'], 'all_detrendRaw_pre', '-v7.3');

% post only
save([ms_save_dir filesep 'all_binary_post.mat' ], 'all_binary_post', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_post.mat'], 'all_RawTraces_post', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_post.mat'], 'all_detrendRaw_post', '-v7.3');


% pre REM only
save([ms_save_dir filesep 'all_binary_pre_REM.mat'], 'all_binary_pre_REM', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_pre_REM.mat'], 'all_RawTraces_pre_REM', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_pre_REM.mat'], 'all_detrendRaw_pre_REM', '-v7.3');

% post SW only
save([ms_save_dir filesep 'all_binary_post_REM.mat' ], 'all_binary_post_REM', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_post_REM.mat'], 'all_RawTraces_post_REM', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_post_REM.mat'], 'all_detrendRaw_post_REM', '-v7.3');

% pre REM only
save([ms_save_dir filesep 'all_binary_pre_SW.mat'], 'all_binary_pre_SW', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_pre_SW.mat'], 'all_RawTraces_pre_SW', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_pre_SW.mat'], 'all_detrendRaw_pre_SW', '-v7.3');

% post SW only
save([ms_save_dir filesep 'all_binary_post_SW.mat' ], 'all_binary_post_SW', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_post_SW.mat'], 'all_RawTraces_post_SW', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_post_SW.mat'], 'all_detrendRaw_post_SW', '-v7.3');


%% save the LFP arrays if needed


if exist('csc', 'var')
    lfp_mat_dir = [ms_save_dir filesep 'LFP_mats'];
    mkdir(lfp_mat_dir)
    f_list = {'d', 't', 'LG','MG','HG', 'UHG', 'Rip'};
    
    for iF = 1:length(f_list)
        % pre
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre.mat'], ['all_' f_list{iF} '_pre'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_freq_pre.mat'], ['all_' f_list{iF} '_freq_pre'], '-v7.3');

        % pre REM
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_REM.mat'], ['all_' f_list{iF} '_pre_REM'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_freq_pre_REM.mat'], ['all_' f_list{iF} '_freq_pre_REM'], '-v7.3');

        % pre SW
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_pre_SW.mat'], ['all_' f_list{iF} '_pre_SW'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_freq_pre_SW.mat'], ['all_' f_list{iF} '_freq_pre_SW'], '-v7.3');

        
        %post
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_post.mat'], ['all_' f_list{iF} '_post'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_freq_post.mat'], ['all_' f_list{iF} '_freq_post'], '-v7.3');
        
        % post REM
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_REM.mat'], ['all_' f_list{iF} '_post_REM'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_freq_post_REM.mat'], ['all_' f_list{iF} '_freq_post_REM'], '-v7.3');
        
        % post SW
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_post_SW.mat'], ['all_' f_list{iF} '_post_SW'], '-v7.3');
        save([lfp_mat_dir filesep 'all_' f_list{iF} '_freq_post_SW.mat'], ['all_' f_list{iF} '_freq_post_SW'], '-v7.3');

    end
    
end
%% save a file with the zscore config
% save([ms_save_dir filesep ms_fname_save '.mat'], 'ms_seg_resize', '-v7.3')
cfgs = [];
cfgs.z_theshold = z_threshold;
cfgs.date = datestr(date, 'yyyy_mm_dd');
if exist('csc', 'var')
    cfgs.filters.d = cfg_d;
    cfgs.filters.t = cfg_t;
    cfgs.filters.LG = cfg_lg;
    cfgs.filters.MG = cfg_mg;
    cfgs.filters.HG = cfg_hg;
    cfgs.filters.UHG = cfg_uhg;
    cfgs.filters.Rip = cfg_rip;
end
save([ms_save_dir filesep 'cfgs_z_' strrep(num2str(z_threshold), '.','p') '_on_' cfgs.date   '.mat'], 'cfgs', '-v7.3')

close all
end % end function

 %     if length(ms_seg.time) ~= length(ms_seg.NLX_evt.t{end})
    %         warning(['Segment # ' num2str(iSeg) ', length of time and NLX_events do not match. Skipping...'])
    %         continue
    %     end
    % check for inactive cells and remove from ms.SFPs just using sum of
    % binary > 0;
    %     %% if using the NLX_csc then filter here.  Checks if the csc input var is a struct (as a real csc would be). if not use the built in NLX_csc
    %     if exist('csc','var') && ~isstruct(csc)
    %
    %         if ~isfield(ms_seg, 'NLX_csc')
    %             error('ms file does not have an ''NLX_csc'' field')
    %         end
    %         this_csc = ms_seg.NLX_csc;
    %         LFP_chan = contains(this_csc.label, 'LFP'); % find the lfp channel.
    %
    %
    %
    %                         % add in the amplitude for each miniscope timestamp
    %         % get the indicies for each miniscope
    %         these_idx = nearest_idx3(ms_seg.NLX_evt.t{end},this_csc.tvec');
    %
    %         % remove other channels
    %         this_csc.data(~LFP_chan,:) = [];
    %         this_csc.label = this_csc.label{LFP_chan};
    %         this_csc.cfg.hdr = this_csc.cfg.hdr(LFP_chan);
    %
    %         cfg_d = [];
    %         % filters
    %         cfg_d.type = 'fdesign'; %Cheby1 is sharper than butter
    %         cfg_d.f  = [1 5]; % broad, could use 150-200?
    %         cfg_d.order = 8; %type filter order (fine for this f range)
    %         cfg_d.display_filter = 0; % use this to see the fvtool
    %
    %         delta = FilterLFP(cfg_d, this_csc);
    %
    %
    %         cfg_t = [];
    %         % filters
    %         cfg_t.type = 'cheby1'; %Cheby1 is sharper than butter
    %         cfg_t.f  = [6 11]; % broad, could use 150-200?
    %         cfg_t.order = 3; %type filter order (fine for this f range)
    %         cfg_t.display_filter = 0; % use this to see the fvtool
    %
    %         theta = FilterLFP(cfg_t, this_csc);
    %
    %
    %         cfg_lg = [];
    %         cfg_lg.check = 0; % plot checks.
    %         % filters
    %         cfg_lg.type = 'cheby1'; %Cheby1 is sharper than butter
    %         cfg_lg.f  = [30 55]; %
    %         cfg_lg.order = 5; %type filter order (fine for this f range)
    %         cfg_lg.display_filter =0; % use this to see the fvtool
    %
    %         LG = FilterLFP(cfg_lg, this_csc);
    %
    %         cfg_hg = [];
    %         % filters
    %         cfg_hg.type = 'cheby1'; %Cheby1 is sharper than butter
    %         cfg_hg.f  = [70 90]; %
    %         cfg_hg.order = 5; %type filter order (fine for this f range)
    %         cfg_hg.display_filter =0; % use this to see the fvtool
    %
    %         HG = FilterLFP(cfg_hg, this_csc);
    %
    %         cfg_rip = [];
    %         % filters
    %         cfg_rip.type = 'butter'; %Cheby1 is sharper than butter
    %         cfg_rip.f  = [120 250]; %
    %         cfg_rip.order = 4; %type filter order (fine for this f range)
    %         cfg_rip.display_filter =0; % use this to see the fvtool
    %
    %         Ripple = FilterLFP(cfg_rip, this_csc);
    %
    %         % get the sampling rate
    %         Fs = this_csc.cfg.hdr{1}.SamplingFrequency;
    %         %     Fs = floor(1/(0.001*(mode(diff((ms_seg_resize.time{1}))))));
    %
    %         D_amp = smooth(abs(hilbert(delta.data)), floor(Fs*0.1));
    %         T_amp = smooth(abs(hilbert(theta.data)), floor(Fs*0.1));
    %         LG_amp = smooth(abs(hilbert(LG.data)), floor(Fs*0.1));
    %         HG_amp = smooth(abs(hilbert(HG.data)), floor(Fs*0.1));
    %         Rip_amp = smooth(abs(hilbert(Ripple.data)), floor(Fs*0.1));
    %
    %
    %
    %         ms_seg.d_amp = D_amp(these_idx);
    %         ms_seg.t_amp = T_amp(these_idx);
    %         ms_seg.LG_amp = LG_amp(these_idx);
    %         ms_seg.HG_amp = HG_amp(these_idx);
    %         ms_seg.Rip_amp = Rip_amp(these_idx);
    %
    %         % check plot
    %         figure(100)
    %         ax(1)=subplot(8,1,1:2);
    %         hold on
    %         c_ord = [linspecer(10);linspecer(10)];
    %         for iC = 1:20
    %             plot(ms_seg.time*0.001, ms_seg.RawTraces(:,iC)+iC, 'color', c_ord(iC,:))
    %         end
    %         xlim([ms_seg.time(1) ms_seg.time(end)]*0.001)
    %
    %         ax(2)= subplot(8,1,3);
    %         hold on
    %
    %         plot(this_csc.tvec- this_csc.tvec(1), this_csc.data,'color', 'k');
    %         plot(this_csc.tvec- this_csc.tvec(1), delta.data,'color', c_ord(1,:));
    %
    %         plot(this_csc.tvec- this_csc.tvec(1), theta.data,'color', c_ord(2,:));
    %
    %         legend({'Raw', 'Delta', 'Theta'})
    %         xlim([0 this_csc.tvec(end)- this_csc.tvec(1)])
    %
    %         ax(3)=subplot(8,1,4);
    %         plot(ms_seg.time*0.001, ms_seg.d_amp,'color', c_ord(1,:));
    %         legend('delta')
    %         xlim([ms_seg.time(1) ms_seg.time(end)]*0.001)
    %
    %         ax(4)=subplot(8,1,5);
    %         plot(ms_seg.time*0.001, ms_seg.t_amp,'color',c_ord(2,:));
    %         legend('theta')
    %         xlim([ms_seg.time(1) ms_seg.time(end)]*0.001)
    %
    %         ax(5)=subplot(8,1,6);
    %         plot(ms_seg.time*0.001, ms_seg.LG_amp,'color',c_ord(3,:));
    %         legend('low gamma')
    %         xlim([ms_seg.time(1) ms_seg.time(end)]*0.001)
    %
    %         ax(6)=subplot(8,1,7);
    %         plot(ms_seg.time*0.001, ms_seg.HG_amp,'color',c_ord(4,:));
    %         legend('high gamma')
    %         xlim([ms_seg.time(1) ms_seg.time(end)]*0.001)
    %
    %         ax(7)=subplot(8,1,8);
    %         plot(ms_seg.time*0.001, ms_seg.Rip_amp,'color',c_ord(5,:));
    %         legend('Ripple')
    %         xlim([ms_seg.time(1) ms_seg.time(end)]*0.001)
    %         linkaxes(ax, 'x')
    %         pause(.5)
    %         close all
    %     elseif  exist('csc','var')