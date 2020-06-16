%% MS_detect_SWR_EV
% loop through subjects/sessions and find lfp events.
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
    PARAMS.data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\Processed place'; % where to find the raw data
    PARAMS.raw_data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LT&sleep'; % raw data location.
    PARAMS.csc_data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LFP'; % where are the LFP files. If blank will look in the same folder as raw_data.
    PARAMS.inter_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\Stats'; % where to put the statistical output .txt
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

%% loop through subjects/sessions.
Subjects = MS_list_dir_names_any(PARAMS.data_dir, {'5'}); % get the subject names.  limited to dir that start with '5' since subjects are 535, 537, and 540
%Subjects{iSub} = Subjects{1} % example.
for iSub = 1:length(Subjects)
    % set the dir for this subject
    this_sub_dir = [PARAMS.data_dir filesep Subjects{iSub}];
    
    % move to dir and get all sessions.  Can specify a specific string to
    % find in the dir. ex: MS_list_dir_names(this_sub_dir, 'LTD')
    cd(this_sub_dir)
    sess_list = MS_list_dir_names_any(this_sub_dir, {Subjects{iSub}}); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    %     sess_list = MS_list_dir_names(this_sub_dir); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    
    
    for iSess = 1:length(sess_list)
        %         iSess = sess_list{4}% picking the 4th session from the list,.
        
        % get the day id (EVA data struct)
        if contains(sess_list{iSess}, 'base') %baseline need to have the base in the title while task sessions do not.
            parts = strsplit(sess_list{iSess}, 'day');
            day_id = ['day' parts{end}(1:6)];
        else
            parts = strsplit(sess_list{iSess}, 'day');
            day_id = ['day' parts{end}(1)];
        end
        
        % crazy line to convert between time formats.
        ms_date = datestr(datenum(strrep(sess_list{iSess}(1:(strfind(sess_list{iSess}, '2019')+3)),'_', '/'), 'MM/dd/yyyy'), 'yyyy-MM-dd');
        
        % find the right CSC dir.  use MS_list_dir_names for specific matches with
        % multiple strings, or MS_list_dir_names_any for any dir names that contain
        % one of the search strings.
        if contains(sess_list{iSess}, 'base')
            [csc_dir_found, csc_dir_fold] = MS_list_dir_names([PARAMS.csc_data_dir filesep day_id(1:4)], {ms_date, day_id});
        else
            [csc_dir_found, csc_dir_fold] = MS_list_dir_names_any([PARAMS.csc_data_dir filesep day_id(1:4)], ms_date);
        end
        
        if ~isempty(csc_dir_found) && length(csc_dir_found) ==1
            fprintf('<strong>%s</strong>: found csc folder at: <strong>%s</strong>\n','MS_segment_data_workflow_sandbox', csc_dir_found{1});
            csc_dir{1} = [csc_dir_fold{1} filesep csc_dir_found{1}];
        elseif ~isempty(csc_dir_found) && length(csc_dir_found) >1
            fprintf('<strong>%s</strong>: Too many CSC folders found\n','MS_segment_data_workflow_sandbox');
            for ii = 1:length(csc_dir_found)
                fprintf('<strong>%s</strong>: found csc folder at: <strong>%s</strong>\n','MS_segment_data_workflow_sandbox', csc_dir_found{ii});
                csc_dir{ii} = [csc_dir_fold{ii} filesep csc_dir_found{ii}];
            end
        else
            error('No csc folder found.');
        end
        
        %         ms_resize_dir = [PARAMS.inter_dir filesep iSub filesep iSess];  %just save the ms_resize struct back into the same place as the ms.mat file.
        %         mkdir(ms_resize_dir);
        
                cd(csc_dir{1})

        
        %% Segment and select the data
        cfg_load = [];
        % for loading the csc
        if strcmpi(Subjects{iSub}, '540')
            cfg_load.fc = {'CSC6.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
            %         elseif strcmpi(Subjects{iSub}, 'PV1069')
            %             cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
            %         elseif strcmpi(Subjects{iSub}, 'PV1043')
            %             cfg_seg.csc.fc = {'CSC1.ncs','CSC6.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        elseif strcmpi(Subjects{iSub}, '537')
            cfg_load.fc = {'CSC7.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        elseif strcmpi(Subjects{iSub}, '535') % csc5 has great SWD and Theta but poor SWR.
            cfg_load.fc = {'CSC7.ncs'}; %
        else
            cfg_load.fc = {'CSC6.ncs'}; % use csc files from Keys if you have them. Alternatively, just use the actual names as: {'CSC1.ncs', 'CSC5.ncs'};
        end
        cfg_load.label = {'LFP'}; % custom naming for each channel.
        cfg_load.desired_sampling_frequency = 2000;
        
        
        cfg.check.saveas = '.png'; % if this is specified it will save the figures in this format. Comment to bypass saving.
        
        %% load data
        csc = MS_LoadCSC(cfg_load);
        
        %% Spike Wave Discharges (SWDs)
        
        cfg_swd = [];
        cfg_swd.check = 0; % plot checks.
        % filters
        cfg_swd.filt.type = 'cheby1'; %Cheby1 is sharper than butter
        cfg_swd.filt.f  = [220 800]; % based on EV suggestion
        cfg_swd.filt.order = 4; %type filter order (fine for this f range)
        cfg_swd.filt.display_filter = 0; % use this to see the fvtool
        
        % use kernel
        %                 cfg_swd.kernel.samples = 60;
        %                 cfg_swd.kernel.sd = 20;
        
        
        % artifact removal
        %                 cfg_swd.artif_det = [];  % toggle artifact removal.
        %                 cfg_swd.artif_det.method = 'zscore';
        %                 cfg_swd.artif_det.threshold = 5;
        %                 cfg_swd.artif_det.dcn = '>';
        % cfg_swd.artif_det.minlen = 0.01;
        
        % detection
        cfg_swd.threshold = 6;% in sd
        cfg_swd.method = 'zscore';
        cfg_swd.min_len = 0.005;
        cfg_swd.merge_thr = 0.01;
        % restrictions
        cfg_swd.max_len = [];
        %                 cfg_swd.max_len.operation = '<';
        %                 cfg_swd.max_len.threshold = .1;
        %                 cfg_swd.nCycles = 3; % number of cycles
        
        % variaence
        %                 cfg_swd.var = [];
        %                 cfg_swd.var.operation = '>';
        %                 cfg_swd.var.threshold = 0.05;
        
        SWD_evts = MS_get_LFP_events_sandbox(cfg_swd, csc);
        
        close all
        cfg_plot.display = 'tsd';
        PlotTSDfromIV(cfg_plot, SWD_evts, csc)
        %                 saveas(gcf, 'SWD_evts.fig');
        % Keep SWD indices for removal in other channels
        cfg_re = [];
        cfg_re.d = [-0.1 0.1];
        SWD_evts = ResizeIV(cfg_re, SWD_evts);
        SWD_idx = TSD_getidx2(csc,SWD_evts); % if error, try TSD_getidx (slower)
        
        %% run event detection for SWRs
        
        cfg_swr = [];
        cfg_swr.check = 0; % plot checks.
        cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
        cfg_swr.filt.f  = [120 250]; % broad, could use 150-200?
        cfg_swr.filt.order = 4; %type filter order (fine for this f range)
        cfg_swr.filt.display_filter = 0; % use this to see the fvtool
        
        % artifact removal (for SWDs that got away)
        %                 cfg_swr.artif_det.threshold = 2.5;
        %                 cfg_swr.artif_det.method = 'zscore';
        %                 cfg_swr.artif_det.rm_len = 0.25;
        %                 cfg_swr.artif_det.dcn = '>';
        
        % smoothing
        cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
        cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;
        
        % detection
        cfg_swr.threshold = 2.5;% in sd
        cfg_swr.method = 'zscore';
        cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
        cfg_swr.merge_thr = 0.01; %merge events that are within 20ms of each other.
        cfg_swr.nan_idx = SWD_idx; % where are any nans, say from excluding artifacts, other events...
        
        % restrictions
        cfg_swr.max_len = [];
        cfg_swr.max_len.operation = '<';
        cfg_swr.max_len.threshold = .1;
        
        %                 cfg_swr.min_len = [];
        %                 cfg_swr.min_len.operation = '<';
        %                 cfg_swr.min_len.threshold = .2;
        cfg_swr.nCycles = 20; % number of cycles
        cfg_swr.nCycles_operation = '=<'; % number of cycles
        
        % variaence
        cfg_swr.var = [];
        cfg_swr.var.operation = '<';
        cfg_swr.var.threshold = 1;
        
        SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc);
        
        cfg_plot.display = 'tsd';
        %                 PlotTSDfromIV(cfg_plot, SWR_evts, csc)
        %                 pause(2)
        close all
        
        
        %% get theta blocks (mostly to exclude false-positive SWR
        cfg_th = [];
        cfg_th.check = 1; % plot checks.
        % filters
        cfg_th.filt.type = 'cheby1'; %Cheby1 is sharper than butter
        cfg_th.filt.f  = [6 9]; % broad, could use 150-200?
        cfg_th.filt.order = 3; %type filter order (fine for this f range)
        cfg_th.filt.display_filter = 0; % use this to see the fvtool
        
        % use kernel
        cfg_th.kernel.samples =csc.cfg.hdr{1}.SamplingFrequency/100;
        cfg_th.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;
        
        % artifact removal
        %                 cfg_th.artif_det = [];  % toggle artifact removal.
        %                 cfg_th.artif_det.threshold = 2.5;
        %                 cfg_th.artif_det.method = 'zscore';
        %                 cfg_th.artif_det.rm_len = 0.25;
        %                 cfg_th.artif_det.dcn = '>';
        %                 cfg_swd.artif_det.minlen = 0.01;
        
        % detection
        cfg_th.threshold = .75;% in sd
        cfg_th.method = 'zscore';
        cfg_th.min_len = 1;
        cfg_th.merge_thr = 0.05;
        cfg_th.nan_idx = SWD_idx; % where are any nans, say from excluding artifacts, other events...
        
        % restrictions
        %                 cfg_th.max_len = [];
        %                 cfg_th.max_len.operation = '<';
        %                 cfg_th.max_len.threshold = .1;
        %                 cfg_th.nCycles = 2; % number of cycles
        
        % variaence
        %                 cfg_th.var = [];
        %                 cfg_th.var.operation = '<';
        %                 cfg_th.var.threshold = 1;
        
        th_evts = MS_get_LFP_events_sandbox(cfg_th, csc);
        
        
        %% remove SWRs that overlap with SWD;
        SWR_evts = DifferenceIV([], SWR_evts, SWD_evts);
        %update the SWR that do not occur during theta events
        SWR_evts = DifferenceIV([], SWR_evts, th_evts);
        
        cfg_plot = [];
        cfg_plot.display = 'iv';
        cfg_plot.mode = 'center';
        cfg_plot.width = 0.2;
        cfg_plot.target = csc.label{1};
        cfg_plot.title = 'var_raw';
        PlotTSDfromIV(cfg_plot,SWR_evts,csc);
        
        figure(1)
        saveas(gcf, 'SWR_evts_examples', 'png');
        close all
        
        %% gamma events (low)
        
        cfg_lg = [];
        cfg_lg.check = 0; % plot checks.
        % filters
        cfg_lg.filt.type = 'cheby1'; %Cheby1 is sharper than butter
        cfg_lg.filt.f  = [30 90]; % broad, could use 150-200?
        cfg_lg.filt.order = 5; %type filter order (fine for this f range)
        cfg_lg.filt.display_filter = 0; % use this to see the fvtool
        
        % use kernel
        %                 cfg_lg.kernel.samples = 20;
        %                 cfg_lg.kernel.sd = 20;
        
        % artifact removal
        cfg_lg.artif_det = [];  % toggle artifact removal.
        cfg_lg.artif_det.threshold = 2.5;
        cfg_lg.artif_det.method = 'zscore';
        cfg_lg.artif_det.rm_len = 0.25;
        cfg_lg.artif_det.dcn = '>';
        %                 cfg_swd.artif_det.minlen = 0.01;
        
        % detection
        cfg_lg.threshold = 1;% in sd
        cfg_lg.method = 'zscore';
        cfg_lg.min_len = 0.015;
        cfg_lg.merge_thr = 0.025;
        % restrictions
        cfg_lg.max_len = [];
        cfg_lg.max_len.operation = '<';
        cfg_lg.max_len.threshold = .1;
        cfg_lg.nCycles = 3; % number of cycles
        
        % variaence
        cfg_lg.var = [];
        cfg_lg.var.operation = '<';
        cfg_lg.var.threshold = 1;
        
        lg_evts = MS_get_LFP_events_sandbox(cfg_lg, csc);
        %                 ms_seg_resize.low_gamma_evts{iB} = DifferenceIV([],lg_events, ms_seg_resize.SWD_evts{iB}) ;
        close all
        
        %%
        str_len = length([Subjects{iSub},sess_list{iSess}]);
        fprintf('<strong>%s %s :</strong>',Subjects{iSub},sess_list{iSess})
        fprintf(repmat(' ', 1,abs(str_len-23)))
        fprintf('SWD = %3d %3.2f/min | SWR = %3d %3.2f/min | theta = %3d %3.2f/min | low_G = %3d %3.2f/min \n',...
            length(SWD_evts.tstart), length(SWD_evts.tstart)/(csc.tvec(end)-csc.tvec(1)),...
            length(SWR_evts.tstart),length(SWR_evts.tstart)/(csc.tvec(end)-csc.tvec(1)),...
            length(th_evts.tstart),length(th_evts.tstart)/(csc.tvec(end)-csc.tvec(1)),...
            length(lg_evts.tstart),length(lg_evts.tstart)/(csc.tvec(end)-csc.tvec(1)))
        
        %% collect the events and cfgs and convert time to idx.
        
        events = [];
        events.SWD.iv = SWD_evts;
        events.SWD.cfg = cfg_swd;
        events.SWD.center_idx =  nearest_idx3(IVcenters(SWD_evts), csc.tvec);
        events.SWD.tstart_idx =  nearest_idx3(SWD_evts.tstart, csc.tvec);
        events.SWD.tend_idx =  nearest_idx3(SWD_evts.tend, csc.tvec);
        
        events.SWR.iv = SWD_evts;
        events.SWR.cfg = cfg_swd;
        events.SWR.center_idx =  nearest_idx3(IVcenters(SWR_evts), csc.tvec);
        events.SWR.tstart_idx =  nearest_idx3(SWR_evts.tstart, csc.tvec);
        events.SWR.tend_idx =  nearest_idx3(SWR_evts.tend, csc.tvec);
        
        save('SWR_evts.mat', 'events', '-v7.3')
        
        clear csc csc_dir
    end % session
end % subject