%% sandbox_EV_SWD_SWR_examples
inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\EV_inter\'; % where to save everything.
data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\final_analysis';

ft_dir = 'C:\Users\ecarm\Documents\GitHub\fieldtrip';

Subs = MS_list_dir_names(data_dir);


Fs = 1500; %


%% collect all the ripples and SWD LFPs for each animal
if exist([inter_dir 'SWR_SWD_LFPs.mat'], 'file')
    disp('File found, loading...')
    load([inter_dir 'SWR_SWD_LFPs.mat']);
else
    fprintf('No file found. Collecting data...\n')
    for iSub = length(Subs):-1:1
        sessions = MS_list_dir_names([data_dir filesep Subs{iSub}]);
        
        for iSess = length(sessions):-1:1
            % get the Fs from a csc file.
            session_names{iSub, iSess} = sessions{iSess};
            
            cd([data_dir filesep Subs{iSub} filesep sessions{iSess}]);
            fprintf('Collecting data from %s...\n', cd)
            
            
            ripple = load('ripple_1.mat');
            ripple_1.LFP{iSub, iSess} = ripple.ripple_lfp;
            %         ripple_1.data{iSub, iSess} = ripple_data;
            ripple_1.ind{iSub, iSess} = ripple.ripple_ind;
            
            clear ripple
            
            % get the ripple 2
            ripple = load('ripple_2.mat');
            ripple_2.LFP{iSub, iSess} = ripple.ripple_lfp;
            %         ripple_2.data{iSub, iSess} = ripple_data;
            ripple_2.ind{iSub, iSess} = ripple.ripple_ind;
            
            clear ripple
            
            % get SWD 1
            swd = load('swd_1.mat');
            swd_1.LFP{iSub, iSess} = swd.swd_lfp;
            %         swd_1.data{iSub, iSess} = swd_data;
            swd_1.ind{iSub, iSess} = swd.swd_ind;
            
            clear swd
            
            % get the swd 2
            swd = load('swd_2.mat');
            swd_2.LFP{iSub, iSess} = swd.swd_lfp;
            %         swd_2.data{iSub, iSess} = swd_data;
            swd_2.ind{iSub, iSess} = swd.swd_ind;
            
            clear swd
        end
    end
    % for tracking sessions
    session_table = cell2table(session_names, 'VariableNames', {'Day1' 'Day3' 'Day5'}, 'rownames', Subs);
    ripple_1.session_table = session_table;
    ripple_2.session_table = session_table;
    swd_1.session_table = session_table;
    swd_2.session_table = session_table;
    
    save([inter_dir 'SWR_SWD_LFPs.mat'], 'ripple_1', 'ripple_2', 'swd_1', 'swd_2');
end
%%  compile LFP events across subjects.
clear all*
subjects = ripple_1.session_table.Properties.RowNames;
for iSub = 1:length(subjects)
    % empty for each subject.
    
    
    %     all_ripple_1.(['S_' subjects{iSub}]) = [];
    %     all_ripple_2.(['S_' subjects{iSub}]) = [];
    %     all_swd_1.(['S_' subjects{iSub}]) = [];
    %     all_swd_2.(['S_' subjects{iSub}]) = [];
    
    for iSess = 1:size(ripple_1.session_table,2)
        fprintf('collecting: %s\n', ripple_1.session_table{iSub, iSess}{1});
        all_ripple_1{iSub, iSess} = [];
        all_ripple_2{iSub, iSess} = [];
        all_swd_1{iSub, iSess} = [];
        all_swd_2{iSub, iSess} = [];
        
        % ripples 1
        for iB = 1:length(ripple_1.LFP{iSub,iSess})
            if isempty(ripple_1.LFP{iSub,iSess}{iB})
                continue
            else
                all_ripple_1{iSub, iSess} = [all_ripple_1{iSub, iSess}, ripple_1.LFP{iSub,iSess}{iB}];
            end
        end
        
        % ripples 2
        for iB = 1:length(ripple_2.LFP{iSub,iSess})
            if isempty(ripple_2.LFP{iSub,iSess}{iB})
                continue
            else
                all_ripple_2{iSub, iSess} = [all_ripple_2{iSub, iSess}, ripple_2.LFP{iSub,iSess}{iB}];
            end
        end
        
        
        % SWD 1
        for iB = 1:length(swd_1.LFP{iSub,iSess})
            if isempty(swd_1.LFP{iSub,iSess}{iB})
                continue
            else
                all_swd_1{iSub, iSess} = [all_swd_1{iSub, iSess}, swd_1.LFP{iSub,iSess}{iB}];
            end
        end
        
        % SWD 2
        for iB = 1:length(swd_2.LFP{iSub,iSess})
            if isempty(swd_2.LFP{iSub,iSess}{iB})
                continue
            else
                all_swd_2{iSub, iSess} = [all_swd_2{iSub, iSess}, swd_2.LFP{iSub,iSess}{iB}];
            end
        end
        
    end
end


%% use FieldTrip to compute PETSpec



addpath(ft_dir)
ft_defaults


% load
load([inter_dir 'ft_hdr']);

this_data = [all_ripple_1{:}];
fig_name = 'SWR_pre';
% split cell into trials
for ii = size(this_data,2):-1:1
    these_trls{ii} = this_data(:,ii)';
    these_times{ii} = -1:1/Fs:1;
end

% build a pseudo trl file
data_trl = [];
data_trl.label = {'CSC'};
data_trl.fsample = Fs;
data_trl.cfg = [];
data_trl.hdr = hdr;
data_trl.trial = these_trls;
data_trl.time = these_times;

% use FT to compute the Spec
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = data_trl.label{1};
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';

if contains(lower(fig_name), 'swd')
    cfg.foi          = 150:2:750; % frequencies of interest
    cfg.toi          = -.05:1/data_trl.fsample:0.05; % times of interest
else
    cfg.foi          = 100:2:200; % frequencies of interest
    cfg.toi          = -.2:1/data_trl.fsample:0.2; % times of interest
end
cfg.t_ftimwin    = ones(size(cfg.foi)).*0.05;%20./cfg.foi;  % window size: fixed at 0.5s

cfg.pad          = 'nextpow2'; % recommened by FT to make FFT more efficient.
cfg.feedback     = 'no'; % might supress verbose output, possible speed improvement.
TFR = ft_freqanalysis(cfg, data_trl);

xrange= cfg.toi; 
% track config for plotting.
freq_params_str = sprintf('Spec using %0.0d swrs. Method: %s, Taper: %s', length(trl),cfg.method, cfg.taper);


%%
%         figure(1)
%         subplot(2,3,[1 2 4 5])
cfg = [];
cfg.channel      = data_trl.label{1};
cfg.baseline     = [xrange(1) -.03];
cfg.baselinetype = 'relative';
cfg.title = freq_params_str;
ft_singleplotTFR(cfg, TFR);
xlabel('time (s)');
ylabel('frequency (Hz)')
% add the mean LFP trace.
cut_idx = nearest_idx([xrange(1) xrange(end)], data_trl.time{1});
this_win = this_data(cut_idx(1):cut_idx(2),:);
hold on
if contains(lower(fig_name), 'swd')
    xfac = 40;
else
    xfac = 60;
end
offset = round(median(TFR.freq));
plot(data_trl.time{1}(cut_idx(1):cut_idx(2)), (nanmean(this_win,2).*xfac)+offset, 'w', 'linewidth', 1)
plot(data_trl.time{1}(cut_idx(1):cut_idx(2)), ((nanmean(this_win,2) +std(this_win,[],2)).*xfac)+offset, '--','color', [.8 .8 .8 .8], 'linewidth', .5)
plot(data_trl.time{1}(cut_idx(1):cut_idx(2)), ((nanmean(this_win,2) -std(this_win,[],2)).*xfac)+offset, '--','color', [.8 .8 .8 .8], 'linewidth', .5)



saveas(gcf, [inter_dir filesep 'Spec_' fig_name], 'fig')
saveas(gcf, [inter_dir filesep 'Spec_' fig_name], 'png')
saveas(gcf, [inter_dir filesep 'Spec_' fig_name], 'eps')




