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



subjects = {'pv1043', 'pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254'};
sessions = {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDSwitch'}; 
% 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1069\HATD1'

LFP_dir = 'K:\Jisoo_Project\LFP data\Jisoo';

decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Decoding'; 

ft_dir = 'C:\Users\ecarm\Documents\GitHub\fieldtrip'; 
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
REM_table = array2table(round(all_REM, 2), 'VariableNames', sessions, 'RowNames', subjects);
SWS_table = array2table(round(all_SWS, 2), 'VariableNames', sessions, 'RowNames', subjects);


%% generate a replay event triggered spectrogram for each session
addpath(ft_dir)
ft_defaults;

for iF = 1:length(f_names)
    
    fprintf('<strong>%s</strong>: loading ms_resize from <strong>%s</strong>\n', mfilename, f_names{iF});
    warning off; load([f_names{iF} filesep 'ms_resize.mat']); warning on;
    
    % get the session and subject IDs
    parts = strsplit(f_names{iF},  filesep);
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
    
    % load the CSC files
     % find the right csc folder
    cd(LFP_dir)
    
    this_LFP_dir = MS_list_dir_names(cd, {subject, csc_type});
    
    cd(this_LFP_dir{1});
    
    % get the pre/post times
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
    

    % get the channels to load from the pre-cut ms_seg data.
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
    % restrict the data to the post recording only to avoid FT issues with
    % discontinuous. 
    csc = restrict(csc, EVT.t{S_rec_idx}(post_S_rec_idx), csc.tvec(end)); 
    
    % load the postREM binary for checking lengths of compiled post REM
    % tvec and post REM binary mat
    load([f_names{iF} filesep 'all_binary_post_REM.mat']); 
    
    % get the replay events
    frame = load([decode_dir filesep subject filesep type filesep '1000shuffling.mat'], 'Final_start_frame');
    
    all_post_REM_time = [];
    for ii = 1:length(ms_seg_resize.file_names)
        if strcmpi(ms_seg_resize.hypno_label{ii}, 'REM') && strcmpi(ms_seg_resize.pre_post{ii}, 'Post') 
            if length(ms_seg_resize.NLX_evt{ii}.t{end}) ~= size(ms_seg_resize.time{ii},1)
                this_time = linspace(ms_seg_resize.NLX_evt{ii}.t{end}(1), ms_seg_resize.NLX_evt{ii}.t{end}(end), size(ms_seg_resize.time{ii},1)); 
            else
                this_time = ms_seg_resize.NLX_evt{ii}.t{end}; 
            end
            all_post_REM_time = [all_post_REM_time; this_time'];
        end
    end
    
    if length(all_post_REM_time) ~= size(all_binary_post_REM, 1)
        error('compiled post REM tvec (%.0f) and saved binary (%.0f) differ in length',length(all_post_REM_time),size(all_binary_post_REM, 1) )
        
    end
    
    R_times = all_post_REM_time(frame.Final_start_frame); 
    
      figure(find(contains(sessions, session))*100 + find(contains(sessions, session)))
    Triggered_Spec_FT(csc, R_times, [subject '-' type '  n = ' num2str(length(R_times))])
    hh = hline([6 12], '--r');
    hh(1).Color = [0.91 0.28 0.28]; 
    hh(2).Color = [0.91 0.28 0.28]; 
    xlim([-2 2])
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    c = colorbar;
    caxis([-8 8])
    ylabel(c, 'zscore'); 
    SetFigure([], gcf)
    set(gcf, 'position', [680 416 750 550])


    
end


%% 
for iF = 1:length(ms_seg_resize.file_names)
   fprintf('%s: duration = %.1fsec\n',  ms_seg_resize.file_names{iF}, (ms_seg_resize.time{iF}(end) - ms_seg_resize.time{iF}(1))/1000)
end

%% generate some videos.
% raw_dir = 'K:\Jisoo_Project\RawData\pv1252\11_18_2021_pv1252_HATD1'; 

raw_dir = 'K:\Jisoo_Project\RawData\pv1254\11_13_2021_pv1254_LTD1'; 
ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1254\LTD1';
decoding_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1254\LTD1';
out_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Decoding_data\Videos'; 
% MS_JC_example_REM_video(raw_dir, ms_dir, decoding_dir ,out_dir, [4581,4586,4936,5336])

% MS_JC_example_REM_video(raw_dir, ms_dir, decoding_dir ,out_dir,
% [7791,7796]); % 1252 Hatswitch

MS_JC_example_REM_video(raw_dir, ms_dir, decoding_dir ,out_dir, [1981,2271,2621])


%% try center of mass calculation

% get replays
for iE = length(Final_start_frame):-1:1
events(:,iE) = decoding.REM_decoded_position(Final_start_frame(iE):Final_start_frame(iE)+14); 
end

[N,X] = hist(events(:,iE), length(decoding.bin_centers_vector)); 

x = regionprops(true(size(events(:,iE))),events(:,iE),  'WeightedCentroid');


matrix=events(:,iE)/sum(events(:,iE));
[m,n]=size(matrix);
[I,J]=ndgrid(1:m,1:n);
   centroid=[dot(I(:),matrix(:)),  dot(J(:),matrix(:))]
