%% MS_segment_data_workflow_sandbox
%
%   This workflow has been designed to allow users to specify a data
%   directory for neuralynx recordings ('csc_dir') and miniscope data
%   (ms_dir) which will be used to:
%
%       1) load the miniscope and neuralynx data
%               key functions: MS_collect_timestamps, MS_get_hypno_label,
%               MS_Remove_trace, MS_LoadCSC
%  
%       2) run initial segmentation and alignment based on the TTLs send from the miniscope DAQ to 
%       the neuralynx system as timestamps in the Events.nev file in the csc_dir.
%               key functions:  MS_extract_NLX_blocks_sandbox,
%               MS_segment_ms_sandbox, MS_remove_data_sandbox, restrict, MS_append_data_sandbox
%
%       3) Filter the neuralynx csc channels (one emg channel, N other LFP
%       channels (would only recommend using one good quality channel. For
%       a quick inspection of all the csc files as a PSD run MS_Quick_psd
%       in the csc_dir). Currently uses delta (2-5Hz), theta (6-11Hz), 
%       and the delta./theta & theta./emg ratios for later visualization. 
%               key functions: FilterLFP (vdmLab)
%
%       4) Visualize the calcium traces and LFP/EMG for each of the
%       recording segments. MS_plot_ca_nlx will show the Ca trace along
%       with the emg/lfp channels.  MS_plot_spec_resize will plot a
%       spectrogram of the raw LFP as well as the raw LFP & EMG, and the
%       filtered bands and their ratios. The user will then be prompted to
%       select a start and an end cutoff using the ginput command (first
%       click on figure for the start and second for the end cutoffs. If
%       these values seem good, then hit 'enter' and enter, OR hit 'backspace' 
%       to flag for removal, OR hit anything else to  select the values again.  
%       If the user selects outside of the start and end of the plot then the 
%       cutoff values will be NaN and skipped in MS_resize_segments. Finally 
%       the data can be saved as an 'ms_resize' structure which contains the 
%       segmented and resized ms data with the nlx data appended to the struct (also segmented).
%               key functions: MS_plot_ca_nlx, MS_select_cutoffs, MS_resize_segments
%
%
%       General information:
%           - This is still in the development phase so lots of variables
%           were chosen based on visual inspections.  If you see something
%           you would like implemented add it as an issue on 
%           https://github.com/ecarmichael/CEH2/issues
%           
%           - I would encourage every user to contribute.  If you want to
%           make edits without changing any of my stuff simply make your
%           own branch on git. If you think something you have written
%           would be a good addition or if you fix a bug in here then just
%           make a pull request. I am happy to get feedback and
%           collaboration will make this code better.
%
%           - I use a common framework when I make functions.  I make use
%           of configuration structures ('cfg') on many functions which
%           contain default parameters that can be overwritten by users 
%           when calling the function.   
%
%               example:  [overwrite the plot_type parameter to use '3d' plotting] 
%               cfg_check = []; % good practice to clear this before filling it in
%               cfg_check.plot_type = '3d'; % use the '3d' plot instead of the default '2d' for the Ca trace
%
%               MS_plot_ca_nlx(cfg_check, ms_seg, res_csc)
%
%
%               example: [just use the defaults]
%               MS_plot_ca_nlx([], ms_seg, res_csc)
%
%
%   EC 2020-02-18   initial version 
%          
%% initialize the system information and add paths to code bases and data. [FILL THIS IN WITH YOUR OWN PATHS]
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
    
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/'; % raw data location.
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'J:\Williams_Lab\Raw_data\JC\7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.raw_data_dir = 'J:\Williams_Lab\JC_Sleep\'; % raw data location.
    PARAMS.inter_dir = 'J:\Williams_Lab\JC_Sleep_inter\'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\JC_Sleep_inter\Stats\'; % where to put the statistical output .txt
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
if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os

%% Locations of the data to be processed.  
% 
% the miniscope data should already be processed with an 'ms.mat' file.
% All of the individual recordings should be in the same folder as the
% ms.dat. example:
%    ms_dir - 
%         ms.mat
%         H10_M20_S2_REM98s     
%         H11_M36_S21_REM160s   
%         H13_M49_S55_REM30s    
%         H14_M7_S4_REM113s

% the nlx files can be in their own folder or in the same one as the ms data 
%   csc_dir - 
%         'CSC1.ncs'
%         'CSC2.ncs'
%         'CSC6.ncs'
%         'CSC7.ncs'
%         'CSC8.ncs'
%         'CheetahLogFile.txt'
%         'Events.nev'


%   You can loop over multiple sessions at this point if you like but this example is just for one session. 
%   to loop over I would make a sess_list = {'7_12_2019_PV1069_LTD5',
%   7_13_2019_PV1069_LTD5',...} and then loop that as 
%
%
%  

% this_dir = dir(PARAMS.raw_data_dir);
% sess_list = [];
% for iSess = 1:length(this_dir)
%     if strcmp(this_dir(iSess).name, '.') || strcmp(this_dir(iSess).name, '..')
%         continue
%     else
%         sess_list{iSess} = this_dir(iSess).name;
%     end
% end
% 
% sess_list = sess_list(~cellfun('isempty',sess_list));

% date_vals = [];
% for iSess = 1:length(sess_list)
%     parts = strsplit(sess_list{iSess},'_'); 
%     id_idx = strfind(parts, 'PV');
%     id_idx = find(~cellfun('isempty', id_idx));
%     
% %     for iP = 1:3
% %         if iP ==1
% %         date_vals = [date_vals parts{iP} ];
% %         else
% %             date_vals = [date_vals '-' parts{iP} ];
% %         end
% %     end
%             
%     
%     
%     Subjects{iSess} = parts{id_idx};
% end


%  for iSess = sess_list
%   ms_dir = [PARAMS.data_raw filesep iSess];
%   csc_dir = [PARAMS.data_raw filesep iSess];

ms_dir = 'J:\Williams_Lab\JC_Sleep\11_23_2019_PV1060_HATD5';
csc_dir = 'J:\Williams_Lab\JC_Sleep\11_23_2019_PV1060_HATD5\2019-11-23_10-10-12_PV1060_HATD5';


    parts = strsplit(ms_dir,filesep);
    parts = strsplit(parts{end}, '_');
    id_idx = strfind(parts, 'PV');
    id_idx = find(~cellfun('isempty', id_idx)); 
    
    this_subject = parts{id_idx}; 
% ms_dir = [PARAMS.raw_data_dir '7_12_2019_PV1069_LTD5']; 
% csc_dir = [PARAMS.raw_data_dir '7_12_2019_PV1069_LTD5'];

parts = strsplit(ms_dir,  filesep);
ms_resize_dir = [PARAMS.inter_dir parts{end}];  %just save the ms_resize struct back into the same place as the ms.mat file. 
mkdir(ms_resize_dir);
%% run the quick PSD script to pick the best csc channels  (only needs to be run once)
% warning this is slow because of the high sampling rate in the csc files. 
cd(csc_dir)
MS_Quick_psd

%% Segment and select the data
cd(ms_dir)
cfg_seg = [];
% for loading the csc
if strcmp(this_subject, 'PV1060')
    cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
else
    cfg_seg.csc.fc = {'CSC1.ncs','CSC7.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
end
cfg_seg.csc.label = {'EMG', 'LFP'}; % custom naming for each channel.
cfg_seg.csc.desired_sampling_frequency = 2000;

% filters
% delta
cfg_seg.filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_seg.filt_d.f  = [1 5];
cfg_seg.filt_d.order = 8; %type filter order
% theta
cfg_seg.filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_seg.filt_t.f  = [6 11];
cfg_seg.filt_t.order = 3; %type filter order

cfg_seg.TS_nlx_match = 1; % when evoked, this will resize the ms data to better match the NLX data if it is off by 1 sample.  [

cfg.check.Ca_type = 'RawTraces'; % what to use for checking the LFP and Ca Trace
cfg.check.saveas = '.png'; % if this is specified it will save the figures in this format. Comment to bypass saving. 


% for resizing
% cfg_seg.resize.segments = [1,2,3,4]; % which segments to process leave
% out to use all events. 
cfg_seg.resize.resize = 1; % use the gui to select the data for resizing. 
cfg_seg.resize.save_key = 'return'; % which key to use for saving the selected cutoffs
% cfg_seg.resize.redo_key = 'downarrow'; % which key for redoing the current cutoffs.
cfg_seg.resize.remove_key = 'backspace'; % which key to flag the current segment for removal.
cfg_seg.resize.spec.win_s = 2^10; % spectrogram window size.
cfg_seg.resize.spec.onverlap = cfg_seg.resize.spec.win_s / 2; % overlap
cfg_seg.resize.spec.freq = 0.5:0.1:80; % frequency range for spectrogram.
cfg_seg.resize.spec.lfp_chan = 2; % which channel to use for the spectrogram. 

ms_resize = MS_Segment_raw(cfg_seg, csc_dir, ms_dir, ms_resize_dir); 

