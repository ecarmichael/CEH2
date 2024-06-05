% Master_J20_content



if strcmp(computer, 'GLNXA64')
    
    codebase_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/home/williamslab/Documents/Github/CEH2';
    oasis_dir = '/home/williamslab/Documents/Github/OASIS_matlab';
    
    code_dir = '/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '/home/williamslab/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = '/home/williamslab';
    
    
elseif strcmp(computer, 'MACA64')
    
    codebase_dir = '/Users/ecar/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/Users/ecar/Documents/Github/CEH2';
    oasis_dir = '/Users/ecar/Documents/Github/OASIS_matlab';
    
    code_dir = '/Users/ecar/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '//Users/ecar/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = '/Users/ecar/';
    
else
    
    codebase_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
    ca_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
    
    code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';
    
    RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = 'C:\Users\ecarm\';
    
end

restoredefaultpath
c_d = cd;
% cd(oasis_dir)
% addpath(genpath(oasis_dir));
% oasis_setup

addpath(genpath(ca_dir));
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(genpath(code_dir))


cd(c_d)


move_thresh  = 9;
bin_size = [.5];

%% compile the intermediate data with Ephys and Ca
comp_data = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\'; 
CA_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\RAW Calcium\Inter';
NLX_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\RAW Calcium\LFP';

Subs = MS_list_dir_names(CA_dir);

for iSub = length(Subs):-1:1
    sessions = MS_list_dir_names([CA_dir filesep Subs{iSub}], '2019');
    
    for iSess = length(sessions):-1:1
        % get the Fs from a csc file.
        session_names{iSub, iSess} = sessions{iSess};
        
        cd([NLX_dir filesep  sessions{iSess}(end-3:end)]);
        
        this_dir  = MS_list_dir_names([NLX_dir filesep  sessions{iSess}(end-3:end)], {Subs{iSub}, 'base1'});
        
        cd(this_dir{1})
        
        fprintf('Collecting NLX from %s...\n', cd)
        
        evts = MS_LoadEvents();
        nlx_t0 = sort([evts.t{contains(evts.label, 'TTL')}]);
        nlx_t0 = nlx_t0(1); 
        
%         nlx_t0 = evts.t{find(contains(evts.label, 'Starting Recording'))}; 
        
        cd([CA_dir filesep Subs{iSub} filesep sessions{iSess}]);
        fprintf('Collecting csc from %s...\n', cd)
        
        warning off
        load('ms_resize.mat');
        warning on
        % collec the useful fields
        
        %%
        all_Raw_SW = []; all_Raw_REM = [];
        all_Bin_SW = []; all_Bin_REM = [];
        all_tvec_SW = []; all_tvec_REM = [];
        all_csc_SW = ms_seg_resize.NLX_csc{1};
        all_csc_SW.data = [];
        all_csc_SW.tvec = [];
        all_csc_SW.label = [];
        all_csc_SW.label{1} = 'LFP';
        all_csc_SW.cfg.hdr = []; 
        all_csc_SW.cfg.hdr{1} = ms_seg_resize.NLX_csc{1}.cfg.hdr{2}; 
        all_csc_REM = all_csc_SW;
%         all_SWD_SW =  ms_seg_resize.SWD_evts{1}; 
%         all_SWD_SW.tstart = []; all_SWD_SW.tend = [];  all_SWD_SW.usr = []; 
%         all_SWD_REM = all_SWD_SW; 
        
        for ii = 1:length(ms_seg_resize.RawTraces)
            this_ms = [];
            this_ms.time = ms_seg_resize.time{ii};
            this_ms.RawTraces = ms_seg_resize.RawTraces{ii};
            this_ms.numNeurons = ms_seg_resize.numNeurons;
            this_ms =  msExtractBinary_detrendTraces(this_ms);
            
            if strcmpi(ms_seg_resize.hypno_label{ii}, 'SW')
                all_Raw_SW = [all_Raw_SW; ms_seg_resize.RawTraces{ii}]; 
                all_Bin_SW = [all_Bin_SW; this_ms.Binary];
                all_tvec_SW = [all_tvec_SW; ms_seg_resize.time{ii}];
                all_csc_SW.tvec = [all_csc_SW.tvec;  ms_seg_resize.NLX_csc{ii}.tvec - nlx_t0];
                all_csc_SW.data = [all_csc_SW.data  ms_seg_resize.NLX_csc{ii}.data(2,:)];
%                 all_SWD_SW.tstart = cat(1,all_SWD_SW.tstart,ms_seg_resize.SWD_evts{ii}.tstart - nlx_t0); % note that constructor guarantees column vectors
%                 all_SWD_SW.tend = cat(1,all_SWD_SW.tend,ms_seg_resize.SWD_evts{ii}.tend - nlx_t0);

            elseif strcmpi(ms_seg_resize.hypno_label{ii}, 'REM')
                all_Raw_REM = [all_Raw_REM; ms_seg_resize.RawTraces{ii}]; 
                all_Bin_REM = [all_Bin_REM; this_ms.Binary];
                all_tvec_REM = [all_tvec_REM; ms_seg_resize.time{ii}];
                all_csc_REM.tvec = [all_csc_REM.tvec;  ms_seg_resize.NLX_csc{ii}.tvec - nlx_t0];
                all_csc_REM.data = [all_csc_REM.data  ms_seg_resize.NLX_csc{ii}.data(2,:)];
%                 all_SWD_REM.tstart = cat(1,all_SWD_REM.tstart,ms_seg_resize.SWD_evts{ii}.tstart - nlx_t0); % note that constructor guarantees column vectors
%                 all_SWD_REM.tend = cat(1,all_SWD_REM.tend,ms_seg_resize.SWD_evts{ii}.tend-nlx_t0);
            end
        end % sleep episode
        
        %% retry the SWD detection
        
        
                cfg_swd = [];
                cfg_swd.check = 1; % plot checks.
                % filters
                cfg_swd.filt.type = 'cheby1'; %Cheby1 is sharper than butter
                cfg_swd.filt.f  = [240 750]; % based on EV suggestion
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
                cfg_swd.threshold = 3;% in sd
                cfg_swd.method = 'zscore';
                cfg_swd.min_len = 0;
                cfg_swd.merge_thr = 0.01;
                % restrictions
                cfg_swd.max_len = [];
%                 cfg_swd.max_len.operation = '<';
%                 cfg_swd.max_len.threshold = .1;
%                 cfg_swd.nCycles = 3; % number of cycles
                
                % variaence
                cfg_swd.var = [];
                cfg_swd.var.operation = '>';
                cfg_swd.var.threshold = 0.05;
                
                all_SWD_REM = MS_get_LFP_events_sandbox(cfg_swd, all_csc_REM);
                
                cfg_swd.mean_filt.operation = '>'; 
                cfg_swd.mean_filt.threshold = 0.00005; 
                all_SWD_SW = MS_get_LFP_events_sandbox(cfg_swd, all_csc_SW);

                close all
%                 cfg_plot.display = 'tsd';
%                 PlotTSDfromIV(cfg_plot, SWD_evts, all_csc_REM)
        
        all_csc_REM.data

        
    end % session
    
end % subject

% 
%         ms_cooc.sleep = [];
%         ms_cooc.sleep.dirName = ms_seg_resize.dirName; 
%         ms_cooc.sleep.time = ms_seg_resize.time;
%         ms_cooc.sleep.RawTraces = ms_seg_resize.RawTraces;
%         ms_cooc.sleep.FiltTraces = ms_seg_resize.FiltTraces;
%         ms_cooc.sleep.pre_post = ms_seg_resize.pre_post;
%         ms_cooc.sleep.NLX_csc = ms_seg_resize.NLX_csc;
%         ms_cooc.sleep.NLX_evt = ms_seg_resize.NLX_evt;
%         ms_cooc.sleep.hypno_label = ms_seg_resize.hypno_label;
%         ms_cooc.sleep.time_labels = ms_seg_resize.time_labels;
%         ms_cooc.sleep.SWD_evts = ms_seg_resize.SWD_evts;
%         ms_cooc.sleep.SWR_evts = ms_seg_resize.SWR_evts;
%         ms_cooc.sleep.resize = ms_seg_resize.resize;

%%  Run again but for EVV data
data_dir = ('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV');
fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\checks';
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')
cd(data_dir);
j20_list = dir('*data*');

method = 'binary';

J_out = []; JP_out = [];
J20_session = []; J20_novel_idx = []; D3_idx = []; D5_idx = [];
for ii = 1:length(j20_list)
    J20_session{ii} = j20_list(ii).name;
    cd(data_dir)
    
    load(j20_list(ii).name)
    %
    %         J_out{ii} = Pipeline_Asmbly(j20_list(ii).name,bin_size, move_thresh, method);
    %     JP_out{ii} = Pipeline_Asmbly_place(j20_list(ii).name,bin_size, move_thresh, method);
    %
    %     % Summary plots
    %         Pipline_Asmbly_plot(J_out{ii}, [fig_dir filesep method]);
    %         close all
    %     Pipline_Asmbly_plot(JP_out{ii}, [fig_dir filesep method filesep 'place']);
    %     close all
    
    if ~isempty(strfind(j20_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        J20_novel_idx(ii) = 1;
        D3_idx(ii) = 0;
    end
    
    if ~isempty(strfind(j20_list(ii).name, 'D3'))
        D3_idx(ii) = 1;
        J20_novel_idx(ii) = 0;
    end
    
    if ~isempty(strfind(j20_list(ii).name, 'D5'))
        J20_novel_idx(ii) = 0;
        D3_idx(ii) = 0;
    end
    
end
J20_novel_idx = logical(J20_novel_idx);
D3_idx = logical(D3_idx);
D5_idx = logical(~D3_idx & ~J20_novel_idx);


if ~isempty(J_out)
    save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\J_out_' method '.mat'], 'J_out')
end
if ~isempty(JP_out)
    save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\JP_out_' method '.mat'], 'JP_out')
end