%% MS_LFP_event_reactivation_sandbox


%% load a nice session

load('D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\537\12_5_2019_537day1\ms_resize.mat'); 

%%
       for iB = 1:length(ms_seg_resize.RawTraces)
           len = ms_seg_resize.NLX_csc{iB}.tvec(end) - ms_seg_resize.NLX_csc{iB}.tvec(1);
           str_len = length(strcat(ms_seg_resize.file_names{iB},ms_seg_resize.pre_post{iB}, ms_seg_resize.hypno_label{iB})); 
            fprintf('<strong>%s %s - %s:</strong>', ms_seg_resize.file_names{iB},ms_seg_resize.pre_post{iB}, ms_seg_resize.hypno_label{iB})
            fprintf(repmat(' ', 1,abs(str_len-23)))
            fprintf('SWD:%3d  %2.1f/s   SWR:%3d  %2.1f/s   low_G:%3d  %2.1f/s \n',...
                length(ms_seg_resize.SWD_evts{iB}.tstart),length(ms_seg_resize.SWD_evts{iB}.tstart)/len,....
                length(ms_seg_resize.SWR_evts{iB}.tstart),length(ms_seg_resize.SWR_evts{iB}.tstart)/len,...
                length(ms_seg_resize.low_gamma_evts{iB}.tstart),length(ms_seg_resize.low_gamma_evts{iB}.tstart)/len)
       end
        fprintf('_______________________________________________________________________________\n')
%% start with one block. 
iSeg = 1;

    keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
    keep_idx =keep_idx(find((keep_idx ~= iSeg)));
    
    cfg_rem = [];
    ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
    
    ms_seg = MS_de_cell(ms_seg);
    
    % binarize the trace
    
    ms_seg = msExtractBinary_detrendTraces(ms_seg);
    
    % check for inactive cells and remove from ms.SFPs just using sum of
    % binary > 0; 
    
    cfg_SFP = [];
    cfg_SFP.fnc = '=='; 
    cfg_SFP.remove_val = 0; 
    ms_seg = MS_update_SFP(cfg_SFP, ms_seg);
    
    csc = ms_seg.NLX_csc;
    SWRs = ms_seg.SWR_evts;
    SWRs.tstart = SWRs.tstart - csc.tvec(1);
    SWRs.tend = SWRs.tend - csc.tvec(1); 
    nlx_evts = ms_seg.NLX_evt;
    nlx_evts.t{end} = nlx_evts.t{end} - csc.tvec(1); 
        csc.tvec = csc.tvec - csc.tvec(1); 

    %% try some stuff
    
    % little checker. 
    traces = 1:50:size(ms_seg.Binary,2);
    figure(11)
    subplot(2,1,1)
    cfg.target = 'LFP';
    PlotTSDfromIV(cfg, SWRs_wide, csc)
%     plot(csc.tvec - csc.tvec(1), csc.data(2,:))
    xlim([csc.tvec(1) csc.tvec(end)])
    subplot(2,1,2)
    for iT  = traces
       hold on 
    plot(ms_seg.time/1000,    ms_seg.Binary(:, iT)+(iT*0.02))
    end
    xlim([ms_seg.time(1)/1000, ms_seg.time(end)/1000])
   
    
    %% 
    swr_centers = IVcenters(SWRs); % convert to centered events;

    
    % make a wider SWRs IV
    
    SWRs_wide = SWRs;
    SWRs_wide.tend = SWRs.tend + 1;
    
% get the time idx that matches the SWR centers (use this if you just want
% one frame before and one after or something.
swr_ms_idx_centers = nearest_idx3(swr_centers, nlx_evts.t{end});

% alternative:
% get the idx for the start and end of the event.
swr_ms_idx_tstart = nearest_idx3(SWRs.tstart, nlx_evts.t{end});

swr_ms_idx_tend = nearest_idx3(SWRs.tend, nlx_evts.t{end});

% initialize some matricies to store the co-activity.
co_mat = NaN(size(ms_seg.Binary,2),size(ms_seg.Binary,2),length(swr_ms_idx_tstart)); % Make an empty matrix for co-activity
corr_mat = NaN(size(ms_seg.Binary,2),size(ms_seg.Binary,2),length(swr_ms_idx_tstart)); % Make an empty matrix for correlations
Q_mat = NaN(size(ms_seg.Binary,2),length(swr_ms_idx_tstart));

SWR_cat_data = [];
idx_win = [-50 50]; % window (in index values) around the event.

%% quick check for Q mat

cfg = [];
cfg.fc = {'TT7_SS_01_Good.t' };
cfg.getTTnumbers = 0;
S = LoadSpikes(cfg);

Q = MakeQfromS([], S);

%% try a Q with Ca binaries
Ca_TS = MS_Binary2TS(ms_seg);
Ca_TS.usr = [];
Ca_1 = Ca_TS; Ca_1.t = []; Ca_1.label = [];
Ca_1.t{1} = Ca_TS.t{1}; Ca_1.label{1} = Ca_TS.label{1}; 
plot(Ca_1)
hold on
plot(ms_seg.time, ms_seg.Binary(:,1))

Ca_Q = MakeQfromS([],Ca_TS);

Ca_Q_SWR = MakeQfromS2([], Ca_TS, SWRs_wide)

%% get a co-occurance measure for the while block

cfg = [];
cfg.nShuffle = 1000;
cfg.useMask = 1;
cfg.outputFormat = 'vectorU';
CC = CoOccurQ(cfg,Ca_Q_SWR);
    %% Seq NMF again?

    
    
    
    
    
    
    