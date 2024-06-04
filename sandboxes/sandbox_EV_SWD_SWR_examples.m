%% sandbox_EV_SWD_SWR_examples
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\EV_inter\'; % where to save everything.
data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\final_analysis';

ft_dir = 'C:\Users\ecarm\Documents\GitHub\fieldtrip';

Subs = MS_list_dir_names(data_dir);


Fs = 1500; %
ca_Fs = 30; 

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
            swd_1.pop{iSub, iSess} = []; 
            for iB = 1:length(swd.swd_data)
                if ~isempty(swd.swd_data{iB})
                    swd_1.pop{iSub, iSess} = [swd_1.pop{iSub, iSess}; swd.swd_data{iB}.popmat];
                end
            end
            %         swd_1.data{iSub, iSess} = swd_data;
            swd_1.ind{iSub, iSess} = swd.swd_ind;
            
            clear swd
            
            % get the swd 2
            swd = load('swd_2.mat');
            swd_2.LFP{iSub, iSess} = swd.swd_lfp;
            swd_2.pop{iSub, iSess} = [];
            for iB = 1:length(swd.swd_data)
                if ~isempty(swd.swd_data{iB})
                    swd_2.pop{iSub, iSess} = [swd_2.pop{iSub, iSess}; swd.swd_data{iB}.popmat];
                end
            end
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
    
    mkdir(inter_dir(1:end-1)); 
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
        
        Ras_ripple_1{iSub, iSess} = [];
        Ras_ripple_2{iSub, iSess} = [];
        Ras_swd_1{iSub, iSess} = [];
        Ras_swd_2{iSub, iSess} = [];
        
        % ripples 1
        for iB = 1:length(ripple_1.LFP{iSub,iSess})
            if isempty(ripple_1.LFP{iSub,iSess}{iB})
                continue
            else
%                 d_f = diff(ripple_1.LFP{iSub,iSess}{iB}(1,:));
%                 find(d_f == 0)
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
%% generate SWD spec
% load
load([inter_dir 'ft_hdr']);

this_data = []; 
for ii = 1:size(all_swd_2,1)
    for jj = 1:size(all_swd_2,2)
        this_data = [this_data, all_swd_2{ii, jj}];
    end
end
fig_name = 'SWD_post';
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
% data_trl.cfg.trl = trl

% use FT to compute the Spec
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = data_trl.label{1};
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';

if contains(lower(fig_name), 'swd')
    cfg.foi          = 150:5:750; % frequencies of interest
    cfg.toi          = -.5:1/data_trl.fsample:0.5; % times of interest
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
freq_params_str = sprintf('Spec using %0.0d swrs. Method: %s, Taper: %s', length(data_trl.trial),cfg.method, cfg.taper);


%%
% plot the overall popmatrix 
all_pop = []; 
for ii = 1:size(swd_2.pop,1)
   for jj = 1:size(swd_2.pop,2)
    all_pop = [all_pop; swd_2.pop{ii, jj}]; 
   end
end

figure
ca_tvec = ((0:size(all_pop, 2))/ca_Fs) -.5;
imagesc(ca_tvec, 1:size(all_pop,1), all_pop)
set(gca, 'xtick', [-0.5 0 0.5]);
xlabel('time (s)');
ylabel('SWD events')

% plot a PETSpec
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


% move to a subplot; 
figlist=get(groot,'Children');

newfig=figure(200);
tcl=tiledlayout(newfig,'flow');

for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    ax.Parent=tcl;
    ax.Layout.Tile=i;
end

figs2keep = 200;
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));
saveas(gcf, [inter_dir filesep 'Spec_' fig_name], 'fig')
saveas(gcf, [inter_dir filesep 'Spec_' fig_name], 'png')
print(gcf, [inter_dir filesep 'Spec_' fig_name], '-depsc')

%% generate and plot SWR spec and Ca activity

% load
load([inter_dir 'ft_hdr']);

this_data = []; 
for ii = 1:size(all_ripple_2,1)
    for jj = 1:size(all_ripple_2,2)
        this_data = [this_data, all_ripple_2{ii, jj}];
    end
end
fig_name = 'SWR_post';
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
% data_trl.cfg.trl = trl

% use FT to compute the Spec
cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = data_trl.label{1};
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';

if contains(lower(fig_name), 'swd')
    cfg.foi          = 150:5:750; % frequencies of interest
    cfg.toi          = -.5:1/data_trl.fsample:0.5; % times of interest
else
    cfg.foi          = 100:2:200; % frequencies of interest
    cfg.toi          = -.5:1/data_trl.fsample:0.5; % times of interest
end
cfg.t_ftimwin    = ones(size(cfg.foi)).*0.05;%20./cfg.foi;  % window size: fixed at 0.5s

cfg.pad          = 'nextpow2'; % recommened by FT to make FFT more efficient.
cfg.feedback     = 'no'; % might supress verbose output, possible speed improvement.
TFR = ft_freqanalysis(cfg, data_trl);

xrange= cfg.toi; 
% track config for plotting.
if contains(lower(fig_name), 'swd')
freq_params_str = sprintf('Spec using %0.0d SWDs. Method: %s, Taper: %s', length(data_trl.trial),cfg.method, cfg.taper);
else
    freq_params_str = sprintf('Spec using %0.0d SWRs. Method: %s, Taper: %s', length(data_trl.trial),cfg.method, cfg.taper);
end

%%
% plot the overall popmatrix 
all_pop = []; 
for ii = 1:size(ripple_2.pop,1)
   for jj = 1:size(swd_2.pop,2)
    all_pop = [all_pop; swd_2.pop{ii, jj}]; 
   end
end

figure
ca_tvec = ((0:size(all_pop, 2))/ca_Fs) -.5;
imagesc(ca_tvec, 1:size(all_pop,1), all_pop)
set(gca, 'xtick', [-0.5 0 0.5]);
xlabel('time (s)');
ylabel('SWR events')

% plot a PETSpec
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


% move to a subplot; 
figlist=get(groot,'Children');

newfig=figure(200);
tcl=tiledlayout(newfig,'flow');

for i = 1:numel(figlist)
    figure(figlist(i));
    ax=gca;
    ax.Parent=tcl;
    ax.Layout.Tile=i;
end

figs2keep = 200;
all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));
saveas(gcf, [inter_dir filesep 'Spec_' fig_name '_new'], 'fig')
saveas(gcf, [inter_dir filesep 'Spec_' fig_name '_new'], 'png')
print(gcf, [inter_dir filesep 'Spec_' fig_name '_new'], '-depsc')
saveas(gcf, [inter_dir filesep 'Spec_' fig_name '_new'], 'svg')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate a trace image for SWD and SWR

for iSub = length(Subs):-1:1
    sessions = MS_list_dir_names([data_dir filesep Subs{iSub}]);
    
    for iSess = length(sessions):-1:1
        % get the Fs from a csc file.
        session_names{iSub, iSess} = sessions{iSess};
        
        cd([data_dir filesep Subs{iSub} filesep sessions{iSess}]);
        fprintf('Collecting data from %s...\n', cd)
        
        
        ripple = load('ripple_2.mat');
        ripple_1.LFP{iSub, iSess} = ripple.ripple_lfp;
        %         ripple_1.data{iSub, iSess} = ripple_data;
        ripple_1.ind{iSub, iSess} = ripple.ripple_ind;
        
        %             if length(ripple.ripple_lfp) > 2
        %                             for ii = 1:length(ripple.ripple_lfp)
        %
        %                 if size(ripple.ripple_lfp{ii},2) >19
        
        for ii = 1:length(ripple.ripple_lfp)
            for jj = 1:size(ripple.ripple_lfp{ii},2)
                if ii ==3
                    if jj == 20
                        figure(jj)
                        subplot(4,1,1)
                        tvec = ((0:length(ripple.ripple_lfp{ii})-1)/Fs) - 1;
                        plot(tvec,ripple.ripple_lfp{ii}(:,jj));
                        xlim([-0.5 0.5]);
                        
                        subplot(4,1,2:4)
                        ca_tvec = ((0:size(ripple.ripple_data{ii}.mat, 1))/ripple.Fs) -.5;
                        imagesc(ca_tvec, 1:size(ripple.ripple_data{ii}.mat,2), squeeze(ripple.ripple_data{ii}.mat(:,:,jj))');
                        %                 MS_Ca_Raster(squeeze(ripple.ripple_data{ii}.mat(:,:,jj))',ca_tvec,14, [0 0 0])
                        %                 set(gca, 'color', 'w');
                        text(-.5,50, (['rec: ' num2str(ii) ' | evt: ' num2str(jj)]), 'color', 'w');
                        
                        pause()
                        close(jj)
                    end
                end
            end
        end
    end
end
            
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
            
%                         get SWD 2
            swd = load('swd_2.mat');
            swd_2.LFP{iSub, iSess} = swd.swd_lfp;
            swd_2.pop{iSub, iSess} = swd.swd_data;
            swd_2.ind{iSub, iSess} = swd.swd_ind;
            
            for ii = 1:length(swd.swd_lfp)
                for jj = 1:size(swd.swd_lfp{ii},2)
                    figure(jj)
                    subplot(4,1,1)
                    tvec = ((0:length(swd.swd_lfp{ii})-1)/Fs) - 1;
                    plot(tvec,swd.swd_lfp{ii}(:,jj));
                    xlim([-0.5 0.5]);
                    
                    subplot(4,1,2:4)
                    text(-.5,50, (['rec: ' num2str(ii) ' | evt: ' num2str(jj)]));
                    ca_tvec = ((0:size(swd.swd_data{ii}.mat, 1)-1)/swd.Fs) -.5;
                    MS_Ca_Raster(squeeze(swd.swd_data{ii}.mat(:,:,jj))',ca_tvec,14, [0 0 0])
                    set(gca, 'color', 'w');
                    pause(1)
                    close(jj)
                end
            end
            

 
 %% plot all the popmats for SWD
 
 figure
imagesc(swd_1.pop{:})

