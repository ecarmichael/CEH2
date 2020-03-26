%% sandbox for event detection


% set paths
close all
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;

if ismac
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    
elseif strcmp(os, 'GLNXA64')
    
        PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
%     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % raw data location.
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\Inter\PV1060\11_21_2019_PV1060_HATD3'; % where to find the raw data
    PARAMS.raw_data_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo\'; % raw data location.
    PARAMS.csc_data_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo\2019-11-21_09-23-36_PV1060_HATD3'; % where are the LFP files. If blank will look in the same folder as raw_data.
    PARAMS.inter_dir = 'J:\Williams_Lab\evt_temp\'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\evt_temp\'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.raw_data_dir) % move to the data folder

% try the newer NLX loaders for UNIX
[~, d] = version;
if str2double(d(end-3:end)) >2014 && isunix
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os


%% get some LFP data from a full recording? 

% cd(PARAMS.csc_data_dir)
cfg_csc = [];
cfg_csc.fc = {'CSC7.ncs'}; 
cfg_csc.desired_sampling_frequency = 2000; %decimate higher sampling rate data to 2kHz for faster processing and well defined filters. 
csc = MS_LoadCSC(cfg_csc);

%% get a block from a recording? 
% cd(PARAMS.data_dir)
load('ms_resize.mat')
lfp_chan = 2;
% for block_idx = 1:22
 block_idx = 4 ; % sw block
% block_idx = 9; % REM block
csc = ms_seg_resize.NLX_csc{block_idx}; 
csc = MS_TSD_SelectChannel(csc, csc.label{lfp_chan});
% fprintf('Loading ms_seg block <strong>%s</strong>, using csc channel labeled: <strong>%s</strong>...\n', ms_seg_resize.file_names{block_idx}, csc.label{1})
%% run event detection for SWRs

cfg_swr = [];
cfg_swr.check = 1; % plot checks. 
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [140 250]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter = 0; % use this to see the fvtool

% detection
cfg_swr.threshold = 1;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2015
cfg_swr.merge_thr = 0.02; %merge events that are within 20ms of each other.

% restrictions
cfg_swr.max_len = 0.1; 
cfg_swr.nCycles = 3; % number of cycles

SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc); 

cfg_plot.display = 'tsd'; 
PlotTSDfromIV(cfg_plot, SWR_evts, csc)

%% gamma detection
lfp_chan = 2;
block_idx = 9; % REM block
csc = ms_seg_resize.NLX_csc{block_idx}; 
csc = MS_TSD_SelectChannel(csc, csc.label{lfp_chan});

cfg_lg = [];
cfg_lg.check = 1; % plot checks. 
% filters
cfg_lg.filt.type = 'cheby1'; %Cheby1 is sharper than butter
cfg_lg.filt.f  = [40 65]; % broad, could use 150-200?
cfg_lg.filt.order = 5; %type filter order (fine for this f range)
cfg_lg.filt.display_filter = 0; % use this to see the fvtool

% detection
cfg_lg.threshold = 1;% in sd
cfg_lg.method = 'zscore';
cfg_lg.min_len = 0.025;
cfg_lg.merge_thr = 0;
% restrictions
cfg_lg.max_len = []; 
cfg_lg.nCycles = 3; % number of cycles

low_gamma_evts = MS_get_LFP_events_sandbox(cfg_lg, csc); 
% end


%% high gamma detection
cfg_hg = [];
cfg_hg.check = 1; % plot checks. 
% filters
cfg_hg.filt.type = 'cheby1'; %Cheby1 is sharper than butter
cfg_hg.filt.f  = [70 120]; % broad, could use 150-200?
cfg_hg.filt.order = 5; %type filter order (fine for this f range)
cfg_hg.filt.display_filter = 0; % use this to see the fvtool

% detection
cfg_hg.threshold = 1;% in sd
cfg_hg.method = 'zscore';
cfg_hg.min_len = 0.025;
cfg_hg.merge_thr = 0;
% restrictions
cfg_hg.max_len = []; 
cfg_hg.nCycles = 3; % number of cycles

low_gamma_evts = MS_get_LFP_events_sandbox(cfg_hg, csc); 
%% theta detection


cfg_t = [];
cfg_t.check = 1; % plot checks. 
% filters
cfg_t.filt.type = 'cheby1'; %Cheby1 is sharper than butter
cfg_t.filt.f  = [6 12]; % broad, could use 150-200?
cfg_t.filt.order = 3; %type filter order (fine for this f range)
cfg_t.filt.display_filter = 0; % use this to see the fvtool

% detection
cfg_t.threshold = .90;% in sd
cfg_t.method = 'percentile'; 
cfg_t.min_len = .5; % 40ms from Vandecasteele et al. 2015
cfg_t.merge_thr = .5; % merge events that are within 20ms of each other.
% restrictions
cfg_t.max_len = []; 
cfg_t.nCycles = 6; % number of cycles


theta_evts = MS_get_LFP_events_sandbox(cfg_t, csc); 
cfg_plot.display = 'tsd'; 
PlotTSDfromIV(cfg_plot, theta_evts, csc)