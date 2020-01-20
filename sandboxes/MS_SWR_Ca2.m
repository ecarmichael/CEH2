% SWR detection and clustering sandbox

%% add paths

close all
restoredefaultpath
global PARAMS
os = computer; 

if ismac
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.ft_code_dir = '/Users/jericcarmichael/Documents/GitHub/fieldtrip'; % FieldTrip toolbbox (used for spectrogram visualization)

elseif strcmp(os, 'GLNXA64')

    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.chronux_code_dir = '/home/ecarmichael/Documents/chronux/chronux_2_12'; % FieldTrip toolbbox (used for spectrogram visualization)
    PARAMS.ft_code_dir = '/home/ecarmichael/Documents/GitHub/fieldtrip'; % FieldTrip toolbbox (used for spectrogram visualization)

else
 disp('on a PC')
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

% try the newer NLX loaders for UNIX
[~, d] = version;
if str2double(d(end-3:end)) >2014
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os


%% load the Miniscope data

load('ms.mat')

% get the timestamps
raw_data_folder = strsplit(PARAMS.data_dir, filesep);
[TS, TS_name] = MS_collect_timestamps(strjoin([PARAMS.raw_data_dir raw_data_folder(end)], ''));

% compare to TS to ms
fprintf('\n****Comparing TS files to processed miniscope (ms) data\n')
for iT = 1:length(TS)
    if length(TS{iT}.system_clock{1}) == ms.timestamps(iT)
        disp([TS_name{iT}   ':  ' num2str(length(TS{iT}.system_clock{1}))   ' - ms TS: ' num2str(ms.timestamps(iT))  '   ~ ' num2str(length(TS{iT}.system_clock{1}) / TS{iT}.cfg.Fs{1}) 's'])
    else
        warning(['TS do not match ms data' TS{iT}.filename   ':  ' num2str(length(TS{iT}.system_clock{1}))   ' - ms TS: ' num2str(ms.timestamps(iT))])
    end
end

ms_seg = MS_segment_ms_sandbox(ms);

fprintf('\n MS_SWR_Ca2: miniscope data has been segmented into %d individual recording epochs\n method used: %s\n', length(ms_seg.time), ms_seg.format); 

%% Load nlx data

% load the Keys file with all of the experiment details. 
%(can be generated with the 'MS_Write_Keys' function) 
ExpKeys = MS_Load_Keys(); 

% load events
nlx_evts = LoadEvents([]); % get '.nev' file here.  

% load the NLX CSC data (using vandermeer lab code) [todo:replace with own]
cfg_csc = [];
cfg_csc.fc = {'CSC7.ncs'}; % use csc files from Keys. Alternatively, just use the actual names {'CSC1.ncs', 'CSC5.ncs'}; 
% cfg_csc.decimateByFactor = 16;
csc = LoadCSC(cfg_csc); % need to comment out ExpKeys lines in LoadCSC


% extract NLX event epochs
cfg_evt_blocks = [];
cfg_evt_blocks.t_chan = 5;
nlx_evts.t{5} = sort([nlx_evts.t{3}, nlx_evts.t{4}]);
nlx_evts.label{5} = 'merge TTls at 3 and 4';
[evt_blocks, evt_iv, evt_duration] = MS_extract_NLX_blocks_sandbox(cfg_evt_blocks, nlx_evts);


% restrict the detected event blocks using gittering to isolate the fewest
% jumps. 


% compare the TS to the NLX evets

% compare to TS to ms
fprintf('\n****Comparing TS files to processed miniscope (ms) data\n')
if length(TS) ~= length(evt_blocks)
    warning('Number of Timestamp files (%s) does not match the number of detected NLX event blocks (%s)', num2str(length(TS)),num2str(length(evt_blocks)))
end

%% append the NLX data to the ms structure
this_chan = 5;
flag = [];
for iT = 1:length(TS)
    if length(TS{iT}.system_clock{1}) == length(evt_blocks{iT}.t{this_chan})
        disp(['TS' num2str(iT) '-' TS_name{iT} ': ' num2str(length(TS{iT}.system_clock{1}))   'samples, '  num2str(length(TS{iT}.system_clock{1}) / TS{iT}.cfg.Fs{1},3) 'sec at ' num2str(TS{iT}.cfg.Fs{1},3) 'Hz'...
            'NLX: ' num2str(length(evt_blocks{iT}.t{this_chan})) ' samples,' num2str(evt_duration(iT),3) 'at ' num2str(1/mode(diff(evt_blocks{iT}.t{this_chan})),3) 'Hz'])
        res_csc{iT} = restrict(csc, evt_blocks{iT}.t{this_chan}(1), evt_blocks{iT}.t{this_chan}(end));
        res_evt{iT} = restrict(nlx_evts, evt_blocks{iT}.t{this_chan}(1), evt_blocks{iT}.t{this_chan}(end));
        
    else
        warning('TS do not match nlx .nev data. TS# %s  %s samples  - NLX: %s events',...
            num2str(iT), num2str(length(TS{iT}.system_clock{1})), num2str(length(evt_blocks{iT}.t{this_chan})))
        flag = [flag, iT];
        res_csc{iT} = [];
        res_evt{iT} = [];
    end
end
res_csc = res_csc(~cellfun('isempty',res_csc));
res_evt = res_evt(~cellfun('isempty',res_evt));


%% update the ms structure with the NLX data
ms_seg = MS_remove_data_sandbox(ms_seg, [flag]);
fprintf('\n MS_SWR_Ca2: miniscope epoch: %d was flagged for removal\n', flag); 

ms_seg = MS_append_data_sandbox(ms_seg, 'NLX_csc', res_csc, 'NLX_evt', res_evt);
fprintf('\n MS_SWR_Ca2: NLX_csc appended\n');

% clear large variables from workspace for memory. 
clear ms res_csc res_evt flag


%% quick check? 
check = 1; % toggle to skip check plots. 
check_evt = 1;

if check  ==1
    figure(101)
    MS_plot_ca_nlx([], ms, csc)
    
    
    
    
end


%% old block for identifying recording epochs. 
% %% Ca blocks [old: replaced with MS_extract_NLX_blocks_sandbox]
% % identify peaks in  diff(evt.t{3}) marking transitions in the camera TTLs
% 
% t_idx = 3; % which event index to use.  
% 
% peak_threshold =  (mean(diff(nlx_evts.t{t_idx}) +0.05*std(diff(nlx_evts.t{t_idx}))));
% min_dist = 10;
% [Rec_peak, Rec_idx] = findpeaks(diff(nlx_evts.t{t_idx}), 'minpeakheight',peak_threshold, 'minpeakdistance', min_dist);
% fprintf(['\nDetected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(Rec_idx), length(Rec_idx))
% 
% 
% for iRec = 1:length(Rec_idx)
%     if iRec < length(Rec_idx)
%         rec_evt{iRec} = restrict(nlx_evts, nlx_evts.t{t_idx}(Rec_idx(iRec)), nlx_evts.t{t_idx}(Rec_idx(iRec+1))); % restrict the NLX evt struct to ms TTL periods
%     else
%         rec_evt{iRec} = restrict(nlx_evts, nlx_evts.t{t_idx}(Rec_idx(iRec)), nlx_evts.t{t_idx}(end)); % restrict the NLX evt file (last only)
%     end
%     all_rec_evt_len(iRec) = length(rec_evt{iRec}.t{t_idx});
% end
% 
% % find the largest and use that one for now. 
% [~,large_idx] = max(all_rec_evt_len);
% 
% 
% %% use the identified largest recording with what should be MS frame grabs
% csc_res = restrict(csc, rec_evt{large_idx}.t{t_idx}(1), rec_evt{large_idx}.t{t_idx}(end));
% fprintf('\nRestricting to section from events file. Duration: %0.0fsecs = %0.2fmins\n', (rec_evt{large_idx}.t{t_idx}(end) -rec_evt{large_idx}.t{t_idx}(1)), (rec_evt{large_idx}.t{t_idx}(end) -rec_evt{large_idx}.t{t_idx}(1))/60)
% 
% %% use whole data
% 
% csc_res = csc;
% 
% %% initial: use a section that looks like SW sleep [use actual timestamps later but needs MS or TS files]; 
% 
% % temp hack to test dectection
% rec.type = 'ts';
% rec.tstart = 4.413082480877258e+06; % place near the end
% rec.tend = nlx_evts.t{2}(end);
% 
% csc_res = restrict(csc, rec.tstart, rec.tend);
% fprintf('\nRestricting to visually identified section.  Duration: %0.0fsecs = %0.2fmins\n',  rec.tend - rec.tstart,(rec.tend - rec.tstart)/60)

%% basic filtering and thresholding
% mouse SWR parameters are based off of Liu, McAfee, & Heck 2017 https://www.nature.com/articles/s41598-017-09511-8#Sec6
check = 1; % used for visual checks on detected events. 

%set up ripple band 
cfg_filt = [];
cfg_filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_filt.f  = [140 250]; % broad, could use 150-200?
cfg_filt.order = 4; %type filter order (fine for this f range)
cfg_filt.display_filter = 0; % use this to see the fvtool 
csc_ripple = FilterLFP(cfg_filt, csc_res);


% convert to amplitude or power
amp_ripple = csc_ripple; % clone to make things simple and replace


for iChan = 1:size(csc_ripple.data,1)
    amp_ripple.data(iChan,:) = abs(csc_ripple.data(iChan,:));
    % Convolve with a gaussian kernel (improves detection)
    kernel = gausskernel(60,20); % note, units are in samples; for paper Methods, need to specify Gaussian SD in ms
    fprintf('\nGausskernal using 60 samples = %0.0fms with SD = 20 samples (%0.0fms)\n', (60/csc_res.cfg.hdr{1}.SamplingFrequency)*1000, (20/csc_res.cfg.hdr{1}.SamplingFrequency)*1000)
    amp_ripple.data(iChan,:) = conv(amp_ripple.data(iChan,:),kernel,'same');
    amp_ripple.units = 'amplitude';
    % if you want units in power use: 
%     amp_ripple.data(iChan,:) =  amp_ripple.data(iChan,:).^2;
%     amp_ripple.units = 'power';

end
 
if check
    figure(10)
    plot(csc_res.tvec, csc_res.data(1,:),'k',csc_ripple.tvec, csc_ripple.data(1,:), 'r',...
        amp_ripple.tvec, amp_ripple.data(1,:),'b')
    legend({'Raw', '140-250 filt', 'Amp'})
end
    
%% remove large amplitude artifacts before SWR detection


csc_artif = csc_res;
for iChan = 1:size(csc_ripple.data,1)
    csc_artif.data(iChan,:) = abs(csc_artif.data(iChan,:)); % detect artifacts both ways
end

cfg_artif_det = [];
cfg_artif_det.method = 'raw';
cfg_artif_det.threshold = std(csc_artif.data(1,:))*5;
% cfg_artif_det.minlen = 0.01;
cfg_artif_det.target = csc_res.label{1};
evt_artif = TSDtoIV(cfg_artif_det,csc_artif);

cfg_temp = []; cfg_temp.d = [-0.5 0.5];
artif_evts = ResizeIV(cfg_temp,evt_artif);


% plot
if check
    plot(113)
    cfg_plot.display = 'iv'; % tsd, iv
    cfg_plot.target = csc_res.label{1};
    PlotTSDfromIV(cfg_plot,artif_evts,csc_artif);
    hline(cfg_artif_det.threshold )
    pause(3); close all;
end

% zero pad artifacts to improve reliability of subsequent z-scoring
artif_idx = TSD_getidx2(csc_res,evt_artif); % if error, try TSD_getidx (slower)
for iChan = 1:size(csc_ripple.data,1)
    csc_ripple.data(iChan,artif_idx) = 0;
    amp_ripple.data(iChan,artif_idx) = 0; 
end

% plot
if check
%     plot(114)
    hold on
    plot(amp_ripple.tvec, csc_ripple.data(1,:),'g');
    plot(amp_ripple.tvec, amp_ripple.data(1,:),'-k');

    pause(3); close all;
end



fprintf('\n MS_SWR_Ca2: %d large amplitude artifacts detected and zero-padded from csc_ripple.\n',length(artif_evts.tstart));


%% isolate candidate events

% get the thresholds
cfg_detect = [];
cfg_detect.operation = '>';
cfg_detect.dcn = cfg_detect.operation; % b/c odd var naming in TSDtoIV
cfg_detect.method = 'zscore';
cfg_detect.threshold = 2.5;
cfg_detect.target = csc.label{1};
cfg_detect.minlen = 0.020; % 40ms from Vandecasteele et al. 2015
cfg_detect.merge_thr = 0.02; % merge events that are within 20ms of each other. 

[swr_evts,evt_thr] = TSDtoIV(cfg_detect,amp_ripple); 

% % now apply to all data
% cfg_select = [];
% cfg_select.dcn = '>';
% cfg_select.method = 'raw';
% cfg_select.threshold = evt_thr;
% cfg_select.target = 'CSC1.ncs';
% cfg_select.minlen = cfg_detect.minlen;
% 
% [evt_ids,~] = TSDtoIV(cfg_select,amp_ripple);



fprintf('\n MS_SWR_Ca2: %d events detected initially.\n',length(swr_evts.tstart));

if check
    cfg_plot = []; 
    cfg_plot.display = 'iv'; 
    cfg_plot.mode = 'center'; 
    cfg_plot.width = 0.2;
    cfg_plot.target = csc.label{1};

    PlotTSDfromIV(cfg_plot,swr_evts,csc);
    pause(2); close all;
end

    
    
    %% exclude events with insufficient cycles - count how many exist above same threshold as used for detection
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filt;
    swr_evt_out = CountCycles(cfg_cc,csc_res,swr_evts);
    
    % get get the evetns with sufficient cycles. 
    cfg_gc = [];
    cfg_gc.operation = '>=';
    cfg_gc.threshold = 5;
    swr_evt_out = SelectIV(cfg_gc,swr_evt_out,'nCycles');
    fprintf('\n MS_SWR_Ca2: %d events remain after cycle count thresholding (%d cycle minimum).\n',length(swr_evt_out.tstart), cfg_gc.threshold);
    
    %% check for evnts that are too long.
    % add in a user field for the length of the events (currently not used)
    swr_evt_out.usr.evt_len = (swr_evt_out.tend - swr_evt_out.tstart)';
    
    cfg_max_len = [];
    cfg_max_len.operation = '<';
    cfg_max_len.threshold = .1;
    swr_evt_out = SelectIV(cfg_max_len,swr_evt_out,'evt_len');
    
    fprintf('\n MS_SWR_Ca2: %d events remain after event length cutoff (> %d ms removed).\n',length(swr_evt_out.tstart), (cfg_max_len.threshold)*1000);
    
    
    %% check for evnts with high raw varience. 'var_raw' is added as a swr_evt_out.usr field in CountCycles

    cfg_max_len = [];
    cfg_max_len.operation = '<';
    cfg_max_len.threshold = 1;
    swr_evt_out = SelectIV(cfg_max_len,swr_evt_out,'var_raw');
    
    fprintf('\n MS_SWR_Ca2: %d events remain after raw varience thresholding (''var_raw'' > %d removed).\n',length(swr_evt_out.tstart), cfg_max_len.threshold);
    
    %% remove events that cooinside with artifacts.
    swr_evt_out = DifferenceIV([], swr_evt_out, artif_evts);
    
    fprintf('\n MS_SWR_Ca2: %d events remain after removing those co-occuring with artifacts.\n',length(swr_evt_out.tstart));

    %% check again
if check
    cfg_plot = []; 
    cfg_plot.display = 'iv'; 
    cfg_plot.mode = 'center'; 
    cfg_plot.width = 0.2;
    cfg_plot.target = csc.label{1};
    cfg_plot.title = 'var';
    PlotTSDfromIV(cfg_plot,swr_evt_out,csc);
    pause(3); close all;
end

%% make a spectrogram of the average SWR 

% spectrogram method using means. 

% %% Try Chronux?
% addpath(genpath(PARAMS.chronux_code_dir));
% disp('Chronux added to path')

%% convert LFP data in to SWR 'Trials'

% DID NO USE. It was a pain to make this work across platforms and MATLAB
% versions. Great functions but not easy to get mex files to work :/
%% using FieldTrip Toolbox  (https://github.com/fieldtrip) 
% 
addpath(PARAMS.ft_code_dir);

ft_defaults
% 
% fc = {'CSC7.ncs'};
% data = ft_read_neuralynx_interp(fc); used to updae TSDtoFT to give
% correct formating. Works as MS_TSDtoFT. 

% convert data to ft format and turn into trials. 
data_ft = MS_TSDtoFT([], csc_res); % convert to ft format. 

swr_centers = IVcenters(swr_evt_out); % get the center of the swr event. 

cfg_trl = [];
cfg_trl.t = cat(1,swr_centers);
cfg_trl.t = cfg_trl.t - data_ft.hdr.FirstTimeStamp;
cfg_trl.twin = [-1 1];
cfg_trl.hdr = data_ft.hdr;

trl = ft_maketrl(cfg_trl);

cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data_ft); 


cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = data_ft.label{1};
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 20:2:250; % frequencies of interest
cfg.t_ftimwin    = ones(size(cfg.foi)).*0.05;%20./cfg.foi;  % window size: fixed at 0.5s
cfg.toi          = -.2:0.0025:0.2; % times of interest
cfg.pad          = 'nextpow2'; % recommened by FT to make FFT more efficient. 
 
TFR = ft_freqanalysis(cfg, data_trl);

% track config for plotting. 
freq_params_str = sprintf('Spec using %0.0d swrs. Method: %s, Taper: %s', length(trl),cfg.method, cfg.taper);

figure
cfg = []; 
cfg.channel      = data_ft.label{1};
cfg.baseline     = [-1 -.01];
cfg.baselinetype = 'relative';
cfg.title = freq_params_str;
ft_singleplotTFR(cfg, TFR);


%% attempt to use chronux. 
% cfg_trials = [];
% 
% swr_centers = IVcenters(swr_evt_out); % get the center of the swr event. 
% % resize around center. 
% swr_center_iv = iv([swr_centers - 0.05, swr_centers + 0.05]);
% 
% % cfg_trial = []; %cfg_trial.target = csc_res.label{1}; cfg_trial.label = csc_res.label{1}; 
% % swr_trials = AddTSDtoIV(cfg_trial, swr_center_iv, csc_res); 
% 
% 
% % convert data into trials
% for iEvt = length(swr_center_iv.tstart):-1:1
%     
%     data_trials(iEvt,:) = csc_res.data(nearest_idx3(csc_res.data, swr_center_iv.tstart(iEvt)));
%     this_data = restrict(csc_res, swr_center_iv.tstart(iEvt), swr_center_iv.tend(iEvt));
% %     data_trials(iEvt,:) = this_data.data; 
% end
% 
% movingwin=[0.01 0.005]; % set the moving window dimensions
% params.Fs=csc_res.cfg.hdr{1}.SamplingFrequency; % sampling frequency
% params.fpass=[50 300]; % frequencies of interest
% params.tapers=[5 9]; % tapers
% params.trialave=1; % average over trials
% params.err=0; % no error computation
% 
% 
% [S1,t,f] = mtspecgramc(data_trials,movingwin,params); % compute spectrogram
% 
% figure(300)
% plot_matrix(S1,t,f);
% xlabel([]); % plot spectrogram
% caxis([8 28]); colorbar;


