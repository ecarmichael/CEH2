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
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.raw_data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/';
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    
elseif strcmp(os, 'GLNXA64')
    
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    
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
    PARAMS.data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\RawData'; % where to find the raw data
    PARAMS.raw_data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\RawData'; % raw data location.
    PARAMS.csc_data_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo'; % where are the LFP files. If blank will look in the same folder as raw_data.
    PARAMS.inter_dir = 'J:\Williams_Lab\JC_Sleep_inter'; % where to put intermediate files
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
% if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
%     rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
%     addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
%     disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
%     which Nlx2MatCSC
% end

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
Subjects = MS_list_dir_names(PARAMS.raw_data_dir);
for iSub = Subjects
    % set the dir for this subject
    this_sub_dir = [PARAMS.raw_data_dir filesep iSub];
    
    % move to dir and get all sessions.  Can specify a specific string to
    % find in the dir. ex: MS_list_dir_names(this_sub_dir, 'LTD')
    cd(this_sub_dir)
    sess_list = MS_list_dir_names(this_sub_dir, 'PV'); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    %     sess_list = MS_list_dir_names(this_sub_dir); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    
    for iSess = sess_list
        ms_dir = [PARAMS.raw_data_dir filesep iSub filesep iSess];
        
        %         if isempty(PARAMS.csc_data_dir)
        %             csc_dir = ms_dir;
        %         end
        
        % crazy line to convert between time formats.
        ms_date = datestr(datenum(strrep(iSess(1:(strfind(iSess, '20')+5)),'_', '/'), 'MM/dd/yyyy'), 'yyyy-MM-dd');
        
        % TODO add in double check using subject and LTD/HAT
        [csc_dir, csc_dir_fold] = MS_list_dir_names(PARAMS.csc_data_dir, ms_date);
        if ~isempty(csc_dir) && length(csc_dir) ==1
            fprintf('<strong>%s</strong>: found csc folder at: <strong>%s</strong>\n','MS_segment_data_workflow_sandbox', csc_dir{1});
            csc_dir = [csc_dir_fold{1} filesep csc_dir{1}];
        elseif ~isempty(csc_dir) && length(csc_dir) >1
            fprintf('<strong>%s</strong>: Too many CSC folders found\n','MS_segment_data_workflow_sandbox', csc_dir{1});
            for ii = 1:length(csc_dir)
                fprintf('<strong>%s</strong>: found csc folder at: <strong>%s</strong>\n','MS_segment_data_workflow_sandbox', csc_dir{ii});
            end
        else
            error('No csc folder found.','MS_segment_data_workflow_sandbox');
        end
        % find the csc files
        % find the matching folder.
        
        
        
        %         TS_files = dir(fullfile(, '**', '*.ncs'));
        %         csc_dir = TS_files(1).folder;
        
        % ms_dir = 'J:\Williams_Lab\JC_Sleep\11_23_2019_PV1060_HATD5';
        % csc_dir = 'J:\Williams_Lab\JC_Sleep\11_23_2019_PV1060_HATD5\2019-11-23_10-10-12_PV1060_HATD5';
        
        
        ms_resize_dir = [PARAMS.inter_dir filesep iSub filesep iSess];  %just save the ms_resize struct back into the same place as the ms.mat file.
        mkdir(ms_resize_dir);
        %% run the quick PSD script to pick the best csc channels  (only needs to be run once)
        % warning this is slow because of the high sampling rate in the csc files.
        if isempty(dir([csc_dir filesep 'PSD_check*'])) %~exist([csc_dir filesep 'PSD_check.png'], 'file')
            fprintf('<strong>MS_segment_data_workflow_sandbox</strong>: no PSD_check file found.  Running MS_Quick_psd...\n');
            cd(csc_dir)
            MS_Quick_psd([], 'fast')
            pos = get(gcf, 'position');
            set(gcf, 'position', [0, 0, pos(3), pos(4)]);
            pause(5)
            close; clear pos;
        end
        
        %% hardcode dir  JISOO USE THIS! 
        
        
        % 1043
%         ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1043\LTD1'; 
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-06-11_09-02-07_pv1043_LTD1'; 
%         iSess = 'LTD1';
%         iSub = '1043'; 
%         
        % LTD5
%          ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1043\LTD5'; 
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-06-15_10-14-07_PV1043_LTD5'; 
%         iSess = 'LTD1';
%         iSub = '1043'; 
%         
        
%1060
%         ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\LTD1';
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-07-15_09-50-03_PV1060_LTD1';
%         iSess = 'LTD1';

% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\LTD5';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-07-19_10-27-46_PV1060_LTD5';
% iSess = 'LTD5';
%
%         ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATD1';
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-11-19_09-59-43_PV1060_HATD1';
%         iSess = 'HATD1';
%
%         ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATD5';
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-11-23_10-10-12_PV1060_HATD5';
%         iSess = 'HATD5';
%
% 
%         ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATDSwitch';
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-11-26_09-59-11_PV1060_HATSwitch';
%         iSess = 'HATDSwitch';
%         iSub = '1060';
%                 
%                 % 1069
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\LTD1';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-07-08_09-03-55_PV1069_LTD1';
% iSess = 'LTD1'; 

% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\LTD5';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-07-12_09-24-26_PV1069_LTD5';
% iSess = 'LTD5'; 


% % % % % ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATD1';
% % % % % csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-10-14_09-37-25_PV1069_HATD1';
% % % % % iSess = 'HATD1'; 


% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATD5';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-10-18_10-02-44_PV1069_HATD5';
% iSess = 'HATD5'; 


% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATDSwitch';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-10-22_09-35-31_PV1069_HATD6_switch';
% iSess = 'HATD6_switch'; 
% 
%         iSub = 'PV1069'; 
%         ms_resize_dir = ms_dir;


                

          %%%% 1191
%         ms_dir = 'K:\Jisoo_Project\RawData\pv1191\5_19_2021_PV1191_HATD1';
%         csc_dir = 'K:\Jisoo_Project\RawData\LFP data\Jisoo\2021-05-19_09-22-27_pv1191_HATD1';
%         iSess = '5_19_2021_PV1191_HATD1';
%         iSub = 'PV1191';
        
%         ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1191/5_23_2021_PV1191_HATD5';
%         csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-05-23_09-26-10_pv1191_HATD5';
%         iSess = '5_23_2021_PV1191_HATD5';
%         iSub = 'PV1191';

%          %Missing ms.mat file.  Probably needs to be CNMFe'd
%         ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1191/5_25_2021_PV1191_HATDS';
%         csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-05-25_10-11-48_PV1191_HATDSwitch';

         %Missing ms.mat file.  Probably needs to be CNMFe'd
%         ms_dir = 'K:\Jisoo_Project\RawData\pv1191\5_25_2021_PV1191_HATDS';
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2021-05-25_10-11-48_PV1191_HATDSwitch';
%         iSess = '5_25_2021_PV1191_HATDS';
%         iSub = 'PV1191';
        
          %%%% 1192
% %         ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1192/4_17_2021_PV1192_HATD1';
% %         csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-04-17_10-06-00_PV1192_HATD1'; 
%         ms_dir = 'K:\Jisoo_Project\RawData\pv1192\4_17_2021_PV1192_HATD1';
%         csc_dir = 'K:\Jisoo_Project\LFP data\2021-04-17_10-06-00_PV1192_HATD1'; 
%         iSess = '4_17_2021_PV1192_HATD1';
%         iSub = 'PV1192';      
% HATD5
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1192\HATD5';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2021-04-21_09-35-07_PV1192_HATD5';
% iSess = 'HATD5';
% iSub = 'PV1192';


% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1192\HATDSwitch';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2021-04-23_09-46-23_PV1192_HATDSwitch';
% iSess = 'HATDSwitch';
% iSub = 'PV1192';
% 
% ms_resize_dir = ms_dir; 
%         ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1192/4_21_2021_PV1192_HATD5';
%         csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-04-21_09-35-07_PV1192_HATD5';
%         iSess = '4_21_2021_PV1192_HATD5';
%         iSub = 'PV1192';

%         ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1192/4_23_2021_PV1192_HATDSwitch'; 
%         csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-04-23_09-46-23_PV1192_HATDSwitch';
%         iSess = '4_23_2021_PV1192_HATDSwitch';
%         iSub = 'PV1192';


        
%         ms_dir = 'K:\Jisoo_Project\RawData\pv1060\11_19_2019_PV1060_HATD1';
%         csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2019-11-19_09-59-43_PV1060_HATD1';
% %         ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1060/11_19_2019_PV1060_HATD1';
% %         csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2019-11-19_09-59-43_PV1060_HATD1';
%         iSess = '11_19_2019_PV1060_HATD1';
%         iSub = 'PV1060';
        

% 1254
%for LTD1

ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\LTD1';
csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2021-11-13_09-21-27_pv1254_LTD1'; 


iSess = 'LTD1';
iSub = 'PV1254';


% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\LTD5';
% csc_dir = 'K:\Jisoo_Project\LFP data\Jisoo\2021-11-17_09-37-11_pv1254_LTD5'; 
% 
% iSess = 'LTD5' 
% iSub = 'PV1254';

ms_resize_dir = ms_dir;

%for LTD3
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_15_2021_pv1254_LTD3';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-15_09-22-38_pv1254_LTD3'; 
% 
% 
% iSess = '11_15_2021_pv1254_LTD3';
% iSub = 'PV1254';

%for LTD5
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_17_2021_pv1254_LTD5';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-17_09-37-11_pv1254_LTD5'; 
% 
% 
% iSess = '11_17_2021_pv1254_LTD5';
% iSub = 'PV1254';

%for HATD1
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_19_2021_pv1254_HATD1';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-19_09-45-12_pv1254_HATD1'; 
% % 
% % 
% iSess = '11_19_2021_pv1254_HATD1';
% iSub = 'PV1254';

%for HATD3
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_21_2021_pv1254_HATD3';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-21_09-33-39_pv1254_HATD3'; 
% iSess = '11_21_2021_pv1254_HATD3';
% iSub = 'PV1254';

%Error ; 
% Index exceeds the number of array elements (19).
% 
% Error in MS_extract_NLX_blocks_sandbox (line 219)
%                     idx_low =find(diff(temp_evt.t{cfg.t_chan}(1:20)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
% 
% Error in MS_Segment_raw (line 208)
% [evt_blocks, ~, evt_duration] = MS_extract_NLX_blocks_sandbox(cfg.evt, nlx_evts);
% 
% Error in Master_segment_data_workflow (line 544)
%         MS_Segment_raw(cfg_seg, csc_dir, ms_dir, ms_resize_dir);

%for HATD5
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_23_2021_pv1254_HATD5';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-23_09-26-07_pv1254_HATD5'; 
% 
% 
% iSess = '11_23_2021_pv1254_HATD5';
% iSub = 'PV1254';

%for HATD5
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_23_2021_pv1254_HATD5';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-23_09-26-07_pv1254_HATD5'; 
% 
% 
% iSess = '11_23_2021_pv1254_HATD5';
% iSub = 'PV1254';

%for HATDSwitch
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1254/11_25_2021_pv1254_HATDSwitch';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-25_09-21-54_pv1254_HATDSwitch'; 
% 
% 
% iSess = '11_25_2021_pv1254_HATDSwitch';
% iSub = 'PV1254';

%for LTD1- pv1252 - csc1 error?
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_12_2021_pv1252_LTD1';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-12_09-17-50_pv1252_LTD1'; 
% 
% 
% iSess = '11_12_2021_pv1252_LTD1';
% iSub = 'PV1252';


%for LTD5
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_16_2021_pv1252_LTD5';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-16_09-20-53_pv1252_LTD5'; 
% 
% 
% iSess = '11_16_2021_pv1252_LTD5';
% iSub = 'PV1252';

%for HATD1
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_18_2021_pv1252_HATD1';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-18_09-29-54_pv1252_HATD1'; 
% 
% 
% iSess = '11_18_2021_pv1252_HATD1';
% iSub = 'PV1252';

%for HATD5
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_22_2021_pv1252_HATD5';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-22_09-58-41_pv1252_HATD5'; 
% 
% 
% iSess = '11_22_2021_pv1252_HATD5';
% iSub = 'PV1252';

%for HATDSwitch
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_24_2021_pv1252_HATDSwitch';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-24_09-36-57_pv1252_HATDSwitch'; 
% 
% 
% iSess = '11_24_2021_pv1252_HATDSwitch';
% iSub = 'PV1252';

%for LTD3
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_14_2021_pv1252_LTD3';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-14_09-14-49_pv1252_LTD3'; 
% 
% 
% iSess = '11_14_2021_pv1252_LTD3';
% iSub = 'PV1252';

%for HATD3
% ms_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1252/11_20_2021_pv1252_HATD3';
% csc_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/LFP data/Jisoo/2021-11-20_09-31-09_pv1252_HATD3'; 
% 
% 
% iSess = '11_20_2021_pv1252_HATD3';
% iSub = 'PV1252';


% Jisoo, run this section and then run the next section (control + return)

      
        
        
%         ms_resize_dir = [PARAMS.inter_dir filesep iSub filesep iSess];  %just save the ms_resize struct back into the same place as the ms.mat file.
%         mkdir(ms_resize_dir);
        
        %% Segment and select the data
        
        %         % get the session info.
        %         parts = strsplit(ms_dir,filesep);
        %         parts = strsplit(parts{end}, '_');
        % %         id_idx = strfind(parts, 'PV');
        % %         id_idx = find(~cellfun('isempty', id_idx));
        % %
        % %         this_subject = parts{id_idx};
        %
        %         parts = strsplit(ms_dir,  filesep);
        %         this_sess = parts{end};
        
        
        cd(ms_dir)
        cfg_seg = [];
        % for loading the csc
        if strcmpi(iSub, 'PV1060')
            cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        elseif strcmpi(iSub, 'PV1069')
            cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'};
        elseif strcmpi(iSub, 'PV1043')
            cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'};
        elseif strcmpi(iSub, 'PV1191')
            cfg_seg.csc.fc = {'CSC1.ncs','CSC7.ncs'};
        elseif strcmpi(iSub, 'PV1254')
            cfg_seg.csc.fc = {'CSC1.ncs','CSC8.ncs'}; % 7 has higher theta but no gamma/SWR
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
        elseif strcmp(iSess, '10_22_2019_PV1069_HATSwitch')
            cfg_seg.bad_block = [12 16];% This index is based on the number of recording blocks in the NLX evt file.
            cfg_seg.bad_block_name = {'H14_M15_S22_remove', 'H14_M54_S54_remove'};% what is the name of the unwanted recording block folder.
        elseif strcmp(iSess, '6_13_2019_PV1043_LTD3')
            cfg_seg.remove_ts = [];
            cfg_seg.remove_nlx_evt = [1];
        elseif strcmp(iSess, '7_17_2019_PV1060_LTD3')
            cfg_seg.remove_ts = [5];
            cfg_seg.remove_nlx_evt = [5,6];
        elseif strcmp(iSess, '7_10_2019_PV1069_LTD3') %added by Jisoo
            cfg_seg.remove_ts = [];
            cfg_seg.remove_nlx_evt = [1];
        elseif strcmp(iSess, '7_15_2019_PV1060_LTD1') %added by Jisoo
            cfg_seg.remove_ts = [1,2,3,4,5,6,7,8,9]; %Event was not recorded at the beginning of the presleep
            cfg_seg.remove_nlx_evt = [];
        elseif strcmp(iSess, '11_19_2019_PV1060_HATD1') %added by EC.  There is an 'H13_M2_S1' that is long in here.  
            cfg_seg.remove_ts = [9]; %Event was not recorded at the beginning of the presleep
            cfg_seg.remove_nlx_evt = [];
        elseif strcmp(iSess, '5_19_2021_PV1191_HATD1') %added by Jisoo
            cfg_seg.remove_ts = []; %Event was not recorded at the beginning of the presleep
            cfg_seg.remove_nlx_evt = [1,2];
        elseif strcmp(iSess, '11_15_2019_PV1060_Homecage_Sleep') %added by Jisoo
            cfg_seg.remove_ts = [1,2,3,4,5]; %Event was not recorded at the beginning of the presleep
            cfg_seg.remove_nlx_evt = [];
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
        % cfg_seg.resize.segments = [1,2,3,4]; % which segments to process leave
        % out to use all events.
        cfg_seg.resize.resize = 1; % use the gui to select the data for resizing.
        cfg_seg.resize.save_key = 'return'; % which key to use for saving the selected cutoffs
        cfg_seg.resize.unclear_key = 'u'; % which key for redoing the current cutoffs.
        cfg_seg.resize.remove_key = 'backspace'; % which key to flag the current segment for removal.
        cfg_seg.resize.spec.win_s = 2^10; % spectrogram window size.
        cfg_seg.resize.spec.onverlap = cfg_seg.resize.spec.win_s / 2; % overlap
        cfg_seg.resize.spec.freq = 0.5:0.1:80; % frequency range for spectrogram.
        cfg_seg.resize.spec.lfp_chan = 2; % which channel to use for the spectrogram.
        cfg_seg.resize.method = 'wavelet'; 
        %     fprintf('<strong>MS_Segment_raw</strong>: processing session: <strong>%s</strong> ...\n',parts{end});
        
        % run the actual segmentation workflow
%         MS_Segment_raw(cfg_seg, csc_dir, ms_dir, ms_resize_dir);
        
        cfg_csc.fc{1} = cfg_seg.csc.fc{2};
        cfg_csc.desired_sampling_frequency = 2000;
        cd(csc_dir)
        csc = MS_LoadCSC(cfg_csc);
        warning off;
        cd(ms_resize_dir)
        MS_re_binarize_JC(2, ms_resize_dir, ms_resize_dir, 'ms_resize', 'ms_resize', csc);
        % zscore the LFP amp and freq. 
        MS_zscore_LFP_JC('K:\Jisoo_Project\LFP data\Jisoo');
        cd(ms_resize_dir); 
         load('ms_resize.mat'); warning on; 
        [data_out_all, data_out_REM_all, data_out_SWS_all,Threshold, labels] = MS_extract_means_JC([],ms_seg_resize.removed ); %Modified by Jisoo
        % save the within session LFP means. 
        mkdir('AcrossEpisodes');
        Out_all.data_out_all=data_out_all;
        Out_all.data_out_REM_all=data_out_REM_all;
        Out_all.data_out_SWS_all=data_out_SWS_all;
        Out_all.Threshold=Threshold;
        Out_all.labels = labels; 
        save([pwd,'/AcrossEpisodes/Out_all_',num2str(Threshold),'.mat'], 'Out_all')
        
        clearvars -except PARAMS
%         clear csc cfg_csc data_out_all data_out_REM_all data_out_SWS_all Out_all Threshold labels
        %     fprintf('<strong>MS_Segment_raw</strong>: processing session: <strong>%s</strong> complete.\n',parts{end});
        
    end % session
end % subjecti