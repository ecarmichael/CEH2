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

elseif strcmp(os, 'GLNXA64')

    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found

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

%% Load data

% load the Keys file with all of the experiment details. 
%(can be generated with the 'MS_Write_Keys' function) 
ExpKeys = MS_Load_Keys(); 

% load the NLX CSC data (using vandermeer lab code) [todo:replace with own]
cfg_csc = [];
cfg_csc.fc = {ExpKeys.EMG, ExpKeys.goodCSC}; % use csc files from Keys. Alternatively, just use the actual names {'CSC1.ncs', 'CSC5.ncs'}; 
cfg_csc.decimateByFactor = 16;
csc = LoadCSC(cfg_csc); % need to comment out ExpKeys lines in LoadCSC


evt = LoadEvents([]);

%% isolate the two sleep sections with CA imagaing....


%% basic filtering and thresholding
% mouse SWR parameters are based off of Liu, McAfee, & Heck 2017 https://www.nature.com/articles/s41598-017-09511-8#Sec6

%set up ripple band 
cfg_filt = [];
cfg_filt.type = 'cheby1'; %Cheby1 is sharper than butter
cfg_filt.f  = [100 250]; % broad, could use 150-200?
cfg_filt.order = 4; %type filter order (fine for this f range)
cfg_filt.display_filter = 1; % use this to see the fvtool 
csc_ripple = FilterLFP(cfg_filt, csc);


% convert to amplitude or power
amp_ripple = csc_ripple; % clone to make things simple and replace


for iChan = 1:size(csc_ripple.data,1)
    amp_ripple.data(iChan,:) = abs(csc_ripple.data(iChan,:));
    % Convolve with a gaussian kernel (improves detection)
    kernel = gausskernel(60,20); % note, units are in samples; for paper Methods, need to specify Gaussian SD in ms
    fprintf('\nGausskernal using 60 samples = %0.3f ms with SD = 20 samples (%0.3fms)\n', csc.cfg.hdr{1}.SamplingFrequency/60, csc.cfg.hdr{1}.SamplingFrequency/20)
    amp_ripple.data(iChan,:) = conv(amp_ripple.data(iChan,:),kernel,'same');
    amp_ripple.units = 'amplitude';
    % if you want units in power use: 
%     amp_ripple.data(iChan,:) =  amp_ripple.data(iChan,:).^2;
%     amp_ripple.units = 'power';

end
 
if check
    figure(10)
    plot(csc.tvec, csc.data(1,:),'k',csc_ripple.tvec, csc_ripple.data(1,:), 'r',...
        amp_ripple.tvec, amp_ripple.data(1,:),'b')
    legend({'Raw', '150-200 filt', 'Amp'})
end
    
%% remove large amplitude artifacts before SWR detection

check = 1; % used for visual checks on detected events. 

csc_artif = csc;
csc_artif.data(1,:) = abs(csc_artif.data(1,:)); % detect artifacts both ways

cfg_artif_det = [];
cfg_artif_det.method = 'raw';
cfg_artif_det.threshold = std(csc.data(1,:))*4;
cfg_artif_det.minlen = 0;
cfg_artif_det.target = csc.label{1};
evt_artif = TSDtoIV(cfg_artif_det,csc_artif);

cfg_temp = []; cfg_temp.d = [-0.5 0.5];
artif_evts = ResizeIV(cfg_temp,evt_artif);

fprintf('\n MS_SWR_Ca2: %d large amplitude artifacts detected and removed from csc_ripple.\n',length(artif_evts.tstart));

% plot
if check
    cfg_plot.display = 'tsd'; % tsd, iv
    cfg_plot.target = csc.label{1};
    PlotTSDfromIV(cfg_plot,artif_evts,csc);
    pause(2); close all;
end

% zero pad artifacts to improve reliability of subsequent z-scoring
artif_idx = TSD_getidx2(csc,evt_artif); % if error, try TSD_getidx (slower)
for iChan = 1:size(csc_ripple.data,1)
    csc_ripple.data(iChan,artif_idx) = 0;
    amp_ripple.data(iChan,artif_idx) = 0; 
end



%% isolate candidate events

% get the thresholds
cfg_detect = [];
cfg_detect.operation = '>';
cfg_detect.dcn = cfg_detect.operation; % b/c odd var naming in TSDtoIV
cfg_detect.method = 'zscore';
cfg_detect.threshold = 4;
cfg_detect.target = csc.label{1};
cfg_detect.minlen = 0.040; % 40ms from Vandecasteele et al. 2015

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

    
    
    %% exclude events with insufficient gamma cycles - count how many exist above same threshold as used for detection
    cfg_cc = [];
    cfg_cc.threshold_type = 'raw';
    cfg_cc.threshold = evt_thr; % use same threshold as for orignal event detection
    cfg_cc.filter_cfg = cfg_filter;
    evt = CountCycles(cfg_cc,csc,evt);
    
    cfg_cc = [];
    cfg_cc.operation = '>=';
    cfg_cc.threshold = cfg.detect_nCycles-1;
    evt = SelectIV(cfg_cc,evt,'nCycles');
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after cycle count removal.\n',length(evt.tstart),cfg.f_label{iFreq});
 %% check for evnts that are too long. 
    % add in a user field for the length of the events (currently not used)
    evt.usr.evt_len = (evt.tend - evt.tstart)';
   
    cfg_max_len = [];
    cfg_max_len.operation = '<';
    cfg_max_len.threshold = 1;
    evt = SelectIV(cfg_max_len,evt,'evt_len');
    
    fprintf('\n MASTER_CollectGammaEvents: %d %s events remain after max length removal.\n',length(evt.tstart),cfg.f_label{iFreq});

