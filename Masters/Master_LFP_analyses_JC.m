%% Master_LFP_analyses_JC


% get the files to process
f_names  = { 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1043\LTD1', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1043\LTD5',... 
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\LTD1', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\LTD5', ...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATD1', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATD5','C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\LTD1', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\LTD5', ...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATD5','C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1191\HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1192\HATD1', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1192\HATD5','C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1192\HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1252\LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1252\HATD1','C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1252\HATD5', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1252\HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\LTD1', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\HATD1','C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\HATD5', 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\HATDSwitch'};

% 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATD1'

LFP_dir = 'K:\Jisoo_Project\LFP data\Jisoo';


%%  loop over and extract LFP Amp and Freq.

for iF =10%:length(f_names)
    
    % load the ms_seg file to get the LFP filenames and the remove
    ms_resize_dir = f_names{iF};
    
    cd(ms_resize_dir);
    
    % get the session info
    parts = strsplit(cd,  filesep);
    session = parts{end};
    subject = parts{end-1};
    date = parts{end};%(1:10);
    if strcmp(date(end), '_')
        date = date(1:end-1);
    end
    type = strsplit(parts{end}, '_');
    type = type{end};
    if strcmp(type,'HATSwitch') && contains(subject, '1069')
        type = 'HATD6_switch';
    elseif strcmp(type,'HATDSwitch') && contains(subject, '1060')
        type = 'HATSwitch';
    end
    
    
    if exist('LFP_mats_z', 'dir')
        fprintf('<strong>%s</strong>: LFP mats Z found skipping <strong>%s</strong>\n', mfilename, f_names{iF});
        
        log_file.([subject '_' session]) = 'LFP mats z detected skipping';
        continue
    end
    % see if the hypno and ms_seg files exist
    if ~exist([ms_resize_dir filesep 'ms_resize.mat'], 'file')
        warning('No ms_resize.mat found.  Skipping...')
        log_file.([subject '_' session]) = 'No ms_resize.mat found';
        continue
    elseif ~exist([ms_resize_dir filesep 'pREM' filesep 'Hypno.mat'], 'file')
        warning('No Hypno.mat found.  Skipping...');
        log_file.([subject '_' session]) = 'No Hypno.mat found';
        continue
    end
    fprintf('<strong>%s</strong>: loading ms_resize from <strong>%s</strong>\n', mfilename, f_names{iF});
    warning off; load('ms_resize.mat'); warning on;
    
    
    
    % find the right csc folder
    cd(LFP_dir)
    
    this_LFP_dir = MS_list_dir_names(cd, {subject, type});
    
    cd(this_LFP_dir{1});
    
    % get the channels to load from the pre-cut data.
    cfg_load = [];
    for iC = length(ms_seg_resize.NLX_csc{1}.cfg.hdr) % loop over channels.
        cfg_load.fc{iC} = [ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.AcqEntName '.ncs'];
        cfg_load.desired_sampling_frequency = ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.SamplingFrequency;
    end
    cfg_load.fc(1) = [];
    
    % hard code LFP channel
    if strcmp(subject, 'PV1043') && strcmp(type, 'LTD5')
        cfg_load.fc{1} = 'CSC6.ncs';
    end
    
    % load some data.
    fprintf('<strong>%s</strong>: loading csc from <strong>%s</strong>\n', mfilename, cfg_load.fc{1});
    
    csc = MS_LoadCSC(cfg_load);
    
    cd(ms_resize_dir)
    
    MS_re_binarize_JC(2, ms_resize_dir, ms_resize_dir, 'ms_resize', 'ms_resize', csc);
    
    if exist(['pREM' filesep 'Cut_CSC.mat'], 'file')
        fprintf('<strong>%s</strong>: loading CSC_cut <strong>%s</strong>\n', mfilename, cfg_load.fc{1});
        load(['pREM' filesep 'Cut_CSC.mat']);
        csc = CSC_cut;
    end
    
    % zscore the LFP amp and freq.
    MS_zscore_LFP_JC('K:\Jisoo_Project\LFP data\Jisoo', csc);
    
    cd(ms_resize_dir);
    
    [data_out_all, data_out_REM_all, data_out_SWS_all,Threshold, labels] = MS_extract_means_JC([],ms_seg_resize.removed ); %Modified by Jisoo
    
    
    % save the within session LFP means.
    mkdir('AcrossEpisodes');
    Out_all.data_out_all=data_out_all;
    Out_all.data_out_REM_all=data_out_REM_all;
    Out_all.data_out_SWS_all=data_out_SWS_all;
    Out_all.Threshold=Threshold;
    Out_all.labels = labels;
    save([pwd,'/AcrossEpisodes/Out_all_',num2str(Threshold),'.mat'], 'Out_all')
    
    clear date
    log_file.([subject '_' session]) = ['Completed ' date];
    
    
    clearvars -except f_names LFP_dir iF log_file
    close all
end

%% get the REM times; 

subjects = {'pv1043', 'pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254'};
sessions = {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDSwitch'}; 

all_REM = nan(length(subjects), length(sessions)); 
all_SWS = all_REM;
for iF =1:length(f_names)
    
    % load the ms_seg file to get the LFP filenames and the remove
    ms_resize_dir = f_names{iF};
    
    fprintf('<strong>%s</strong>: recomputing REM %% from <strong>%s</strong>\n', mfilename, f_names{iF});
    
    cd(ms_resize_dir);
    
    parts = strsplit(cd,  filesep);
    session = parts{end};
    subject = parts{end-1};
    date = parts{end};%(1:10);
    if strcmp(date(end), '_')
        date = date(1:end-1);
    end
    type = strsplit(parts{end}, '_');
    type = type{end};
    if strcmp(type,'HATDSwitch') && contains(subject, '1060')
        type = 'HATDSwitch';
        csc_type = 'HATSwitch';
    else
        csc_type = type; 
    end
    
    load(['pREM' filesep 'Hypno.mat']);
    
    
    % find the right csc folder and grab the events file
    cd(LFP_dir)
    
    this_LFP_dir = MS_list_dir_names(cd, {subject, csc_type});
    
    cd(this_LFP_dir{1});
    
    EVT = LoadEvents([]);
    S_rec_idx = find(contains(EVT.label, 'Starting Recording')); % get the index for the start recoding cell
    Stop_rec_idx = find(contains(EVT.label, 'Stopping Recording')); % get the index for the start recoding cell
    
    if strcmpi(subject, 'pv1060') && strcmpi(type, 'LTD1')
            pre_S_rec_idx = 2;
        post_S_rec_idx = 3;
    elseif length(EVT.t{S_rec_idx}) ~= 2
        warning('more than two recordings detected.  Fix this later.')
        for iR = length(EVT.t{S_rec_idx}):-1:1
            rec_dur(iR) = EVT.t{Stop_rec_idx}(iR) - EVT.t{S_rec_idx}(iR);
        end
        keep_rec = find((rec_dur/60/60) > 1.5);
        if length(keep_rec) ~= 2
            error('Something is wrong with the CSC. Expected 2 recordings but found %.0f',length(keep_rec));
        end
        pre_S_rec_idx = keep_rec(1);
        post_S_rec_idx = keep_rec(2);       
    else
        pre_S_rec_idx = 1;
        post_S_rec_idx = 2;
    end
    % split for pre-post
    cut_idx = nearest_idx3(EVT.t{S_rec_idx}(post_S_rec_idx),Hypno.tvec);
    Hypno.tvec = Hypno.tvec(cut_idx:end); 
    Hypno.data = Hypno.data(cut_idx:end);

    % get the REM stats. 
    [REM_stats, SWS_stats] = MS_compute_REM_prct_JC(Hypno); 
    
    all_REM(find(contains(subjects, subject)), find(contains(sessions, type))) = REM_stats.percent_sleep; 
    all_SWS(find(contains(subjects, subject)), find(contains(sessions, type))) = SWS_stats.percent_sleep; 

    
end
REM_table = array2table(round(all_REM, 2), 'VariableNames', sessions, 'RowNames', subjects)