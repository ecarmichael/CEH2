% Master_EC_Ca_Maze_preprocess

% code for preprocessing data from calcium mice on the W maze.  Should
% contain an MS data structure and NLX LFP and event data.  

% v1.0   extracts position data from the NLX tracking rather than deep lab
% cut and aligns it with the calcium data. 


%% initialize

close all
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;
if strcmp(os, 'GLNXA64')
    
    if strcmpi(getenv('USERNAME'), 'williamslab')
        PARAMS.raw_data_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData'; % raw data location.
        PARAMS.inter_dir = '/home/williamslab/Dropbox (Williams Lab)/JisooProject2020/2020_Results_aftercutting/Across_episodes/Inter'; % where to put intermediate files
        PARAMS.csc_data_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo'; % where are the LFP files. If blank will look in the same folder as raw_data.
        PARAMS.stats_dir = [PARAMS.inter_dir '/Stats/']; % where to put the statistical output .txt
        PARAMS.code_base_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
        PARAMS.code_CEH2_dir = '/home/williamslab/Documents/Github/CEH2'; % where the multisite repo can be found
        
    elseif strcmpi(getenv('USERNAME'), 'ecarmichael')
        PARAMS.raw_data_dir = '/home/ecarmichael/Dropbox (Williams Lab)'; % raw data location.
        PARAMS.inter_dir     = '/home/ecarmichael/Dropbox (Williams Lab)/JisooProject2020/2020_Results_aftercutting/Across_episodes/Inter'; % where to put intermediate files
        PARAMS.csc_data_dir  = '/mnt/Data'; % where are the LFP files. If blank will look in the same folder as raw_data.
        PARAMS.stats_dir     = [PARAMS.inter_dir '/Stats/']; % where to put the statistical output .txt
        PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
        PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    end
else
    PARAMS.data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\Raw'; % where to find the raw data
    PARAMS.csc_data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\LFP'; % where are the LFP files. If blank will look in the same folder as raw_data.
    PARAMS.inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter'; % where to put intermediate files
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab'; 
end

rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

clear os


%% list the files to process and then loop over each session. 
f_names = dir('*MZD*'); 

for iF = 1:length(f_names)

    data_dir = [f_names(iF).folder filesep f_names(iF).name]; 
    
    csc_dir = [f_names(iF).folder(1:end-4) filesep 'LFP' filesep f_names(iF).name '_LFP']; 
    if ~exist([f_names(iF).folder(1:end-4) filesep 'LFP' filesep f_names(iF).name '_LFP'], 'dir')
        error('CSC dir does not exist');
    end
    % get the subject and session ID
    parts = strsplit(f_names(iF).name, '_');
    iSub = parts{end-1}; 
    iSess = parts{end}; 

    % save the output here. 
    ms_resize_dir = [PARAMS.inter_dir filesep iSub filesep f_names(iF).name]; 
    mkdir(ms_resize_dir); 
    
    
    % set parameters; same as those used by JC
    cfg_seg = [];
    % for loading the csc
    if strcmpi(iSub, 'PV1254')
        cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'}; % 7 has higher theta but no gamma/SWR
    else
        cfg_seg.csc.fc = {'CSC1.ncs','CSC7.ncs'}; % Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
    end
    cfg_seg.csc.label = {'EMG', 'LFP'}; % custom naming for each channel.
    cfg_seg.csc.desired_sampling_frequency = 2000;
    
    % flag known bad recording blocks.  [EC ToDo: make another variable
    % that will remove specific TS or Evt blocks when there is a known
    % discrepency between the MS files and the .nev stamps.
    if strcmp(iSess, '10_18_2019_PV1069_HATD5')
        cfg_seg.bad_block = 19;% This index is based on the number of recording blocks in the NLX evt file.
        cfg_seg.bad_block_name = {'H15_M37_S21_REmove_withoutLFP'};% what is the name of the unwanted recording block folder.
        cfg_seg.remove_ts=[]; % added by Jisoo
        cfg_seg.remove_nlx_evt=[]; % added by Jisoo
    else
        cfg_seg.bad_block = [] ;
        cfg_seg.bad_block_name = {};
        cfg_seg.remove_ts=[]; % added by Jisoo
        cfg_seg.remove_nlx_evt=[]; % added by Jisoo
    end
    
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
    cfg_seg.resize.resize = 1; % use the gui to select the data for resizing.
    cfg_seg.resize.save_key = 'return'; % which key to use for saving the selected cutoffs
    cfg_seg.resize.unclear_key = 'u'; % which key for redoing the current cutoffs.
    cfg_seg.resize.remove_key = 'backspace'; % which key to flag the current segment for removal.
    cfg_seg.resize.spec.win_s = 2^10; % spectrogram window size.
    cfg_seg.resize.spec.onverlap = cfg_seg.resize.spec.win_s / 2; % overlap
    cfg_seg.resize.spec.freq = 0.5:0.1:80; % frequency range for spectrogram.
    cfg_seg.resize.spec.lfp_chan = 2; % which channel to use for the spectrogram.
    cfg_seg.resize.method = 'wavelet';
    
    fprintf('<strong>MS_Segment_raw</strong>: processing session: <strong>%s</strong> ...\n',[iSub '-' iSess]);
    
    %     run the actual segmentation workflow
    addpath(genpath(PARAMS.oasis_dir)); oasis_setup; 
    MS_Segment_raw_EC(cfg_seg, csc_dir, data_dir, ms_resize_dir);
    rmpath(genpath(PARAMS.oasis_dir)); 
    
    cd(ms_resize_dir)
    % extrack the postion data for the
    load('ms_trk.mat');
    
    
    % load the DLC data
    trk_dir = dir([data_dir filesep '*MAZE']);
    if isempty(dir([trk_dir.folder filesep trk_dir.name filesep '*DLC*.csv']))
        error('YOU DIDN"T DLC THIS SESSION YET...!\n')
    end
    [~, behav] = MS_collect_DLC([trk_dir.folder filesep trk_dir.name]);
    
    % convert pixels to cm. 
    behav.position = (behav.position /.666)/10; % convert pixels to cms
    behav.units = 'cm' ;
    behav.t_start =  ms_trk.time(1); % keep the original starting time in case it is needed later.
    behav.time = behav.time + ms_trk.time(1); % correct time to start at zero. 
    if isempty(behav)
        error('behav is empty')
    end
    
    % algin time to matach Ca2+ time 
    if length(ms_trk.time) ~= length(behav.time)
       behav = MS_align_data(behav, ms_trk);
    end
    
    
    % append the idealized coordinates    
    if exist([PARAMS.inter_dir filesep 'Common_CoorD.mat'], 'file') % check for a common coordindate file. 
        load([PARAMS.inter_dir filesep 'Common_CoorD.mat']); % if it exists load it. 
        behav.CoorD_L = Common_CoorD.CoorD_L;
        behav.CoorD_R = Common_CoorD.CoorD_R;
    else
        behav = MS_behav_append_Coord(behav, 'CoorD_L'); % if not then 
        behav = MS_behav_append_Coord(behav, 'CoorD_R'); % if not then
        behav = MS_behav_append_Coord(behav, 'CoorD_M'); % if not then

    end
    
    % append the trial event times from nvt. (assumes NLX zoned tracking
    % was used. If not avilable try using the tracking matlabapp
    % Maze_GUI_sandbox (WIP); 
    load([ms_resize_dir filesep 'Events.mat']); 
    behav = MS_behav_append_MAZE(behav, evt); 
    
    cd(ms_resize_dir)
    save('behav_DLC.mat', 'behav', '-v7.3')


end
