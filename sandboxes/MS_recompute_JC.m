%% MS_recompute_JC
%  Simple loop for Jisoo's CA data.  The core function can be replaced to
%  say rebinarize traces (MS_re_binarize_JC), or load LFP files
%  (MS_LoadCSC)
%
%   EC 2020-07-30   initial version
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
    PARAMS.data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\RawData'; % where to find the raw data
    PARAMS.raw_data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\RawData'; % raw data location.
    PARAMS.csc_data_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo'; % where are the LFP files. If blank will look in the same folder as raw_data.
    PARAMS.inter_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\Inter\Stats'; % where to put the statistical output .txt
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
Subjects = MS_list_dir_names(PARAMS.raw_data_dir);
for iSub = 2:length(Subjects)
    % set the dir for this subject
    this_sub_dir = [PARAMS.raw_data_dir filesep Subjects{iSub}];
    
    % move to dir and get all sessions.  Can specify a specific string to
    % find in the dir. ex: MS_list_dir_names(this_sub_dir, 'LTD')
    cd(this_sub_dir)
%     sess_list = MS_list_dir_names(this_sub_dir, {'LTD', 'HAT'}); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    sess_list = MS_list_dir_names_any(this_sub_dir, {'LTD', 'HAT'}); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'


    for iSess = 1:length(sess_list)
        ms_dir = [PARAMS.raw_data_dir filesep Subjects{iSub} filesep sess_list{iSess}];
        
%         if isempty(PARAMS.csc_data_dir)
%             csc_dir = ms_dir;
%         end
        
        % crazy line to convert between time formats. 
        ms_date = datestr(datenum(strrep(sess_list{iSess}(1:(strfind(sess_list{iSess}, '2019')+3)),'_', '/'), 'MM/dd/yyyy'), 'yyyy-MM-dd');
       
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

        
        ms_inter_dir = [PARAMS.inter_dir filesep Subjects{iSub} filesep sess_list{iSess}];  %just save the ms_resize struct back into the same place as the ms.mat file.
        if ~exist(ms_inter_dir, 'dir') 
            mkdir(ms_inter_dir); % make the inter folder if it does not already exist. 
        end
        %% Select the CSC data to be loaded
        if strcmpi(Subjects{iSub}, 'PV1060')
            cfg_load.fc = {'CSC8.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        elseif strcmpi(Subjects{iSub}, 'PV1069')
            cfg_load.fc = {'CSC7.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        elseif strcmpi(Subjects{iSub}, 'PV1043')
            cfg_load.fc = {'CSC7.ncs'}; % CSC6 was used for resizing given the nice theta. 
        else
            cfg_load.fc = {'CSC7.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        end
        cfg_load.label = { 'LFP'}; % custom naming for each channel.
        cfg_load.desired_sampling_frequency = 2000;        
        
        cfg.check.saveas = '.png'; % if this is specified it will save the figures in this format. Comment to bypass saving.
        
        %% extract and restrict csc blocks 
        
        cd(csc_dir)
        csc = MS_LoadCSC(cfg_load);
        
        cd(ms_inter_dir) % assumes this directory has the pre process/segmented data here.
        if exist([ms_inter_dir filesep 'ms_resize.mat'], 'file')
            MS_re_binarize_JC(2, ms_inter_dir, ms_inter_dir, 'ms_resize', 'ms_resize', csc);
        else
            continue
            warning(sprintf('No ms_resize.mat file can be found in inter dir: <strong>%s</strong>', ms_inter_dir))
        end
        clear csc csc_dir
    end % session
end % subject