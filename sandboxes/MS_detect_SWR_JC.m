%% MS_detect_SWR_JC
% loop through subjects/sessions and find lfp events.  Specific for Jisoo's
% data structure.
%
%
%   EC 2020-06-20   initial version
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
Subjects = MS_list_dir_names(PARAMS.raw_data_dir);
for iSub = length(Subjects)
    % set the dir for this subject
    this_sub_dir = [PARAMS.raw_data_dir filesep Subjects{iSub}];
    
    % move to dir and get all sessions.  Can specify a specific string to
    % find in the dir. ex: MS_list_dir_names(this_sub_dir, 'LTD')
    cd(this_sub_dir)
    sess_list = MS_list_dir_names(this_sub_dir, 'LTD'); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
%     sess_list = MS_list_dir_names(this_sub_dir); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'

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
        mkdir(ms_inter_dir);        
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
        
        %% load data
        
        cd(csc_dir)
        csc = MS_LoadCSC(cfg_load);
        
              %% restrict to section half of sleep (post task)
%         NLX_events = LoadEvents([]);
%         
%         if length(NLX_events.t{1}) ==2
%         csc = restrict(csc, NLX_events.t{1}(2),  NLX_events.t{2}(2));
%         else
%             disp('more than 2 start times')
%         end
        
        %% run event detection for SWRs
        
        cfg_swr = [];
        cfg_swr.check = 0; % plot checks.
        cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
        cfg_swr.filt.f  = [120 250]; % broad, could use 150-200?
        cfg_swr.filt.order = 4; %type filter order (fine for this f range)
        cfg_swr.filt.display_filter =0; % use this to see the fvtool
        
        % artifact removal (for SWDs that got away)
        %                 cfg_swr.artif_det.threshold = 2.5;
        %                 cfg_swr.artif_det.method = 'zscore';
        %                 cfg_swr.artif_det.rm_len = 0.25;
        %                 cfg_swr.artif_det.dcn = '>';
        
        % smoothing
        cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
        cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;
        
        % detection
        cfg_swr.threshold =2;% in sd
        cfg_swr.method = 'zscore';
        cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
        cfg_swr.merge_thr = 0.01; %merge events that are within 20ms of each other.
%         cfg_swr.nan_idx = SWD_idx; % where are any nans, say from excluding artifacts, other events...
        
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
        
%         cfg_plot.display = 'tsd';
%                         PlotTSDfromIV(cfg_plot, SWR_evts, csc)
        %                 pause(2)
        
        g = groot;
        if ~isempty(g.Children)
            figure(1)
            saveas(gcf, 'SWR_evts_examples', 'png');
            close all
        end
        
%         %% get theta blocks (mostly to exclude false-positive SWR
%         cfg_th = [];
%         cfg_th.check = 1; % plot checks.
%         % filters
%         cfg_th.filt.type = 'cheby1'; %Cheby1 is sharper than butter
%         cfg_th.filt.f  = [6 9]; % broad, could use 150-200?
%         cfg_th.filt.order = 3; %type filter order (fine for this f range)
%         cfg_th.filt.display_filter = 0; % use this to see the fvtool
%         
%         % use kernel
%         cfg_th.kernel.samples =csc.cfg.hdr{1}.SamplingFrequency/100;
%         cfg_th.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;
%         
%         % artifact removal
%         %                 cfg_th.artif_det = [];  % toggle artifact removal.
%         %                 cfg_th.artif_det.threshold = 2.5;
%         %                 cfg_th.artif_det.method = 'zscore';
%         %                 cfg_th.artif_det.rm_len = 0.25;
%         %                 cfg_th.artif_det.dcn = '>';
%         %                 cfg_swd.artif_det.minlen = 0.01;
%         
%         % detection
%         cfg_th.threshold = .75;% in sd
%         cfg_th.method = 'zscore';
%         cfg_th.min_len = 1;
%         cfg_th.merge_thr = 0.05;
%         cfg_th.nan_idx = SWD_idx; % where are any nans, say from excluding artifacts, other events...
%         
%         % restrictions
%         %                 cfg_th.max_len = [];
%         %                 cfg_th.max_len.operation = '<';
%         %                 cfg_th.max_len.threshold = .1;
%         %                 cfg_th.nCycles = 2; % number of cycles
%         
%         % variaence
%         %                 cfg_th.var = [];
%         %                 cfg_th.var.operation = '<';
%         %                 cfg_th.var.threshold = 1;
%         
%         th_evts = MS_get_LFP_events_sandbox(cfg_th, csc);
%         
%         
%         figure(1)
%         saveas(gcf, 'SWR_evts_examples', 'png');
%         close all
%         
%         %% gamma events (low)
%         
%         cfg_lg = [];
%         cfg_lg.check = 0; % plot checks.
%         % filters
%         cfg_lg.filt.type = 'cheby1'; %Cheby1 is sharper than butter
%         cfg_lg.filt.f  = [30 90]; % broad, could use 150-200?
%         cfg_lg.filt.order = 5; %type filter order (fine for this f range)
%         cfg_lg.filt.display_filter = 0; % use this to see the fvtool
%         
%         % use kernel
%         %                 cfg_lg.kernel.samples = 20;
%         %                 cfg_lg.kernel.sd = 20;
%         
%         % artifact removal
%         cfg_lg.artif_det = [];  % toggle artifact removal.
%         cfg_lg.artif_det.threshold = 2.5;
%         cfg_lg.artif_det.method = 'zscore';
%         cfg_lg.artif_det.rm_len = 0.25;
%         cfg_lg.artif_det.dcn = '>';
%         %                 cfg_swd.artif_det.minlen = 0.01;
%         
%         % detection
%         cfg_lg.threshold = 1;% in sd
%         cfg_lg.method = 'zscore';
%         cfg_lg.min_len = 0.015;
%         cfg_lg.merge_thr = 0.025;
%         % restrictions
%         cfg_lg.max_len = [];
%         cfg_lg.max_len.operation = '<';
%         cfg_lg.max_len.threshold = .1;
%         cfg_lg.nCycles = 3; % number of cycles
%         
%         % variaence
%         cfg_lg.var = [];
%         cfg_lg.var.operation = '<';
%         cfg_lg.var.threshold = 1;
%         
%         lg_evts = MS_get_LFP_events_sandbox(cfg_lg, csc);
%         %                 ms_seg_resize.low_gamma_evts{iB} = DifferenceIV([],lg_events, ms_seg_resize.SWD_evts{iB}) ;
%         close all
%         
%         %%
%         str_len = length([Subjects{iSub},sess_list{iSess}]);
%         fprintf('<strong>%s %s :</strong>',Subjects{iSub},sess_list{iSess})
%         fprintf(repmat(' ', 1,abs(str_len-23)))
%         fprintf('SWR = %3d %3.2f/min | theta = %3d %3.2f/min | low_G = %3d %3.2f/min \n',...
%             length(SWR_evts.tstart),length(SWR_evts.tstart)/(csc.tvec(end)-csc.tvec(1)),...
%             length(th_evts.tstart),length(th_evts.tstart)/(csc.tvec(end)-csc.tvec(1)),...
%             length(lg_evts.tstart),length(lg_evts.tstart)/(csc.tvec(end)-csc.tvec(1)))
        
        %% collect the events and cfgs and convert time to idx.
        events.SWR.session = sess_list{iSess}; 
        events.SWR.iv = SWR_evts;
        events.SWR.cfg = cfg_swr;
        events.SWR.center_idx =  nearest_idx3(IVcenters(SWR_evts), csc.tvec);
        events.SWR.tstart_idx =  nearest_idx3(SWR_evts.tstart, csc.tvec);
        events.SWR.tend_idx =  nearest_idx3(SWR_evts.tend, csc.tvec);
        
        mkdir([PARAMS.inter_dir filesep 'SWRs']);
        save([PARAMS.inter_dir filesep 'SWRs' filesep sess_list{iSess} '_evts.mat'], 'events', '-v7.3')
        
        %% get the start and stop times for the pre and post parts of the csc. 
        
        NLX_evt = LoadEvents([]);
        
        pre_end = NLX_evt.t{2}(1); 
        post_start = NLX_evt.t{1}(2); 
        
        %% make a cell x time x event array. 
        
              
        % move the miniscope data dir
        cd(ms_inter_dir)
        load('ms_resize.mat');
        
        
        nFrames = 33; % number of frames before and after the center of each event.
        
        % initialize some arrays for later concatenation. 
        all_SWR_activity_pre = [];
        all_SWR_activity_post = [];

%         all_SWR_LFP.data = [];
%         all_SWR_LFP.tvec = []; 
        %% loop through Ca block blocks and get the Ca activity around the SWR.
        
        for iSeg = 1:length(ms_seg_resize.NLX_csc)
            swr_r = restrict(events.SWR.iv, ms_seg_resize.NLX_csc{iSeg}.tvec(1), ms_seg_resize.NLX_csc{iSeg}.tvec(end));
            
            if ~isempty(swr_r.tstart) % if there are no events here, then skip 
                
                swr_centers = IVcenters(swr_r); % get the swr center.
                
                keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
                keep_idx =keep_idx(find((keep_idx ~= iSeg)));
                
                cfg_rem = []; cfg_rem.verbose = 0;
                this_ms = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
                
                this_ms = MS_de_cell(this_ms);
                
                this_ms = MS_msExtractBinary_detrendTraces(this_ms, 2);
                
            
                
                % debugging
                cfg_plot.display = 'tsd';
                cfg_plot.target = 'LFP';
                %                         ca_frames = iv(this_ms.NLX_evt.t{end}-0.001, this_ms.NLX_evt.t{end}+0.001);
                % %                         PlotTSDfromIV(cfg_plot, ca_frames, ms_seg_resize.NLX_csc{iSeg})
                %                         PlotTSDfromIV(cfg_plot, swr_r, ms_seg_resize.NLX_csc{iSeg})
                %                         PlotTSDfromIV(cfg_plot, swr_r, csc)
                
                %
                swr_idx = nearest_idx(swr_centers, this_ms.NLX_evt.t{end}); % get the SWR indicies using NLX events (TTL from frame) which correspond to the ms frame number
                
                
                count = 0; % initialize the counter.
                SWR_activity_pre = [];
                SWR_activity_post = [];

                %                     SWR_LFP_data = []; SWR_LFP_tvec = [];
                
                for iE = length(swr_idx):-1:1
                    % check that the event is not occuring too close to the edge.
                    if swr_idx(iE) < nFrames || swr_idx(iE) > length(this_ms.NLX_evt.t{end})-nFrames
                        continue
                    else
                        count = count+1; % keep a counter for event indexing.  Avoids skipped events.
                        % get and hold the LFP values.
                        %                 this_swr = restrict(csc,this_ms.NLX_evt.t{end}(swr_idx(iE)-nFrames),this_ms.NLX_evt.t{end}(swr_idx(iE)+nFrames));
                        %                 SWR_LFP_data(count,:) = this_swr.data;
                        %                 SWR_LFP_tvec(count,:) = this_swr.tvec;
                        
                        % split into 'pre' or 'post' recording. 
                        if swr_centers(iE) <= pre_end  
                            for iC = this_ms.numNeurons:-1:1
                                SWR_activity_pre(iC, :, count) = this_ms.Binary(swr_idx(iE)-nFrames: swr_idx(iE)+nFrames, iC);
                            end
                        elseif swr_centers(iE) >= post_start
                            for iC = this_ms.numNeurons:-1:1
                                SWR_activity_post(iC, :, count) = this_ms.Binary(swr_idx(iE)-nFrames: swr_idx(iE)+nFrames, iC);
                            end
                        end
                    end
                end
                % append to the master array;
                if swr_centers(end) <= pre_end
                    all_SWR_activity_pre = cat(3,all_SWR_activity_pre, SWR_activity_pre);
                elseif  swr_centers(end) >= post_start
                    all_SWR_activity_post = cat(3,all_SWR_activity_post, SWR_activity_post);
                end
                %             all_SWR_LFP.data = cat(2,all_SWR_LFP.data, SWR_LFP_data);
                
            end % if swr_idx is empty. 
        end  % recording segments/blocks
        %% save the output in the inter folder for this project. 
        
        mkdir([PARAMS.inter_dir filesep 'SWRs']);
        save([PARAMS.inter_dir filesep 'SWRs' filesep sess_list{iSess} '_activity_pre.mat'], 'all_SWR_activity_pre', '-v7.3')
        save([PARAMS.inter_dir filesep 'SWRs' filesep sess_list{iSess} '_activity_post.mat'], 'all_SWR_activity_post', '-v7.3')

        clear csc csc_dir
    end % session
end % subject