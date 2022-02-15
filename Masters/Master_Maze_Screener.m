%% Master_Maze_Screener:  master control script for batch screening cells for spatial informaiton.
%
%
%  To run:  set up your desired directories in the below section and it
%  will process all the sessions in the data_dir. All figures will be saved
%  in the inter_dir along with the intermediate files. Nothing will be
%  written back to the data_dir.
%
%  % this uses PARAMS as a global parameter that can be called across
%  functions.  It is mostly used for tracking directories and colours.
%
% EC 2021-01-02   initial version for subiculum screening.
%
%  TODO:
%   - make the initial data dir selection simple.
%
%% initialize

close all
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;
% parent_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/ck2cre-1359hd/2021_01_30/14_18_06';
% parent_dir = ('/mnt/Data/Behav test and scripts/ck-1361/2021_02_04');
% cd(parent_dir);
% cd('BehavCam_1/')
if strcmp(os, 'GLNXA64')
    
    %% Home
    if strcmpi(getenv('USERNAME'), 'ecarmichael')
        PARAMS.data_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid'; % where to find the raw data
        PARAMS.inter_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/II_inter'; % where to put intermediate files
        PARAMS.stats_dir = '/mnt/Data/Williams_Lab/II_classification/Inter/'; % where to put the statistical output .txt
        PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
        PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
        %% Office
    elseif strcmpi(getenv('USERNAME'), 'williamslab')
        PARAMS.data_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Ingrid (Results_Calcium data 2020)/Project1_NTandCK2comparison'; % where to find the raw data
        PARAMS.inter_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/II_inter'; % where to put intermediate files
        PARAMS.stats_dir = PARAMS.inter_dir; % where to put the statistical output .txt
        PARAMS.code_base_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
        PARAMS.code_CEH2_dir = '/home/williamslab/Documents/Github/CEH2'; % where the multisite repo can be found
    else
        PARAMS.data_dir = '/lustre06/project/6064766/ecar/Maze_ca/inter'; % where to find the raw data
        PARAMS.code_base_dir = '/home/ecar/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
        PARAMS.code_CEH2_dir = '/home/ecar/GitHub/CEH2'; % where the multisite repo can be found
    end    
else
    PARAMS.data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter'; % where to find the rpreprocessed data
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
end

% colours
PARAMS.L_grey = [0.8 0.8 0.8];
PARAMS.D_grey = [0.2 0.2 0.2];
PARAMS.blue = [0.3639    0.5755    0.7484];
PARAMS.red = [0.9153    0.2816    0.2878];
PARAMS.green= [0.4416    0.7490    0.4322];
PARAMS.gold = [1.0000    0.5984    0.2000];

rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

clear d os
beep off % I know when I mess up without that annoying beep, thanks.

% configuration
%general
cfg.method ='Binary';  %'decon'; % can be binary (default) and 'decon'
cfg.binary_thresh = 2; % number of sd for binary thresholding of zscored Ca data.
cfg.split_method = 'time'; % method for splitting session in half.  Can also be 'nTrans' to use number of Ca transients instead.

% place
cfg.p_thres = 0.05; % value for pvalue cut off;
cfg.stability_thres = 0.5; % from van der Veldt 2020
cfg.nShuff = 600;
cfg.p_bin_size = 5 ; % in cm
cfg.split_gaus_sd = 5; % sd for gaussian smoothing of place tuning for split session xcorr.

% speed
cfg.s_bin_size = 1.375;
cfg.s_bins  =  2.5:cfg.s_bin_size:30; % between -2cm/s^2 and 2cm/s^s with 20 bins matches van de Veldt et al. 2020
cfg.s_bins(cfg.s_bins==0) = []; %remove 0 bin.

% acceleration
cfg.accel_bin_size = .1;
cfg.accel_bins  =  -1:cfg.accel_bin_size:1; % between -2cm/s^2 and 2cm/s^s with 20 bins matches van de Veldt et al. 2020
cfg.accel_bins(cfg.accel_bins==0) = []; %remove 0 bin.

% head-direction
cfg.hd_bin_size = 360/15;
cfg.hd_bins = 0:cfg.hd_bin_size:360;

%% navigate the desired directory

% get the list of subjects
sub_list = {};
d = dir;
d=d(~ismember({d.name},{'.','..', '._*'})); % get the folder names and exclude any dir that start with '.'.
for iSess = 1:length(d)
    if ~strcmp(d(iSess).name(1:2), '._') && d(iSess).isdir % exclude any that are autosaves.
        sub_list{end+1} = d(iSess).name; % keep the good folder names.
    end
end

for iSub = 1:length(sub_list) % comment length out and just do iSub = X session.
    
    cd([PARAMS.data_dir filesep sub_list{iSub}])
    
    sess_list = {};
    d = dir;
    d=d(~ismember({d.name},{'.','..', '._*'})); % get the folder names and exclude any dir that start with '.'.
    for iSess = 1:length(d)
        if ~strcmp(d(iSess).name(1:2), '._') && d(iSess).isdir % exclude any that are autosaves.
            sess_list{end+1} = d(iSess).name; % keep the good folder names.
        end
    end
    
    % loop across sessions.
    for iSess = 1:length(sess_list) % loop through sessions for this subject.
        
        cd([PARAMS.data_dir filesep  sub_list{iSub} filesep sess_list{iSess}])
        
        
        % load the trk data
        warning off; load('ms_trk.mat'); warning on;
        
        % get the file info
        parts = strsplit(cd, filesep);
        
        sess_parts = strsplit(strrep(parts{end}, '-', '_'), '_') ;
        
        fname = parts{end};
        f_info = [];
        f_info.subject = parts{end-1};
        f_info.date = datestr([sess_parts{1} '-' sess_parts{2} '-' sess_parts{3}], 'yyyy-mm-dd');
        f_info.task = 'W_MAZE';
        f_info.experimenter = 'EC';
        f_info.fname = fname; % full name.
        time_str = regexp(ms_trk.file_names,'\d*','Match');
        f_info.time = datestr(datestr([time_str{1} ':' time_str{2} ':' time_str{3}]), 'HH:MM:SS');
        
        fprintf('<strong>%s</strong>: loading ms trk from %s %s...\n', mfilename, f_info.fname)
        
        % get the position info from the DLC output file.
        load('behav_DLC.mat')
        
        % preprocess the sessions
%         cfg.cells = [4,68];
        data = MS_compute_place(cfg, ms_trk, behav, f_info);
        
        
        clear ms ms_temp behav
        
        
        % save the output file back to the intermediate dir.
        save([PARAMS.inter_dir filesep strrep([f_info.fname '_' f_info.date '_' strrep(f_info.time, ':','-') '_' strrep(f_info.task, ' ', '_')], '-','_')], 'data', '-v7.3')
        close all
        clear data
    end % end sess
end % end sub


%% collect the processed data and summarize
cd(PARAMS.inter_dir)


%get the subject list
sub_list = {};
d = dir;
d=d(~ismember({d.name},{'.','..', '._*'})); % get the folder names and exclude any dir that start with '.'.
for iSess = 1:length(d)
    if ~strcmp(d(iSess).name(1:2), '._') && d(iSess).isdir % exclude any that are autosaves.
        sub_list{end+1} = d(iSess).name; % keep the good folder names.
    end
end

for iSub = 2:length(sub_list) % comment length out and just do iSub = X session.
    cd(PARAMS.inter_dir)
    
    cd(sub_list{iSub}) % go to the subject folder
    
    sess_list = dir('*.mat');
    all_P_sigs = []; all_P_MI = []; all_sess_P_sigs = NaN(length(sess_list),1); all_sess_P_sigs_split = NaN(length(sess_list),1);  all_P_xcor = [];
    all_S_sigs = []; all_S_MI = []; all_sess_S_sigs = NaN(length(sess_list),1);
    all_fname = cell(length(sess_list),1);
    
    for iSess =  1: length(sess_list)
        
        load(sess_list(iSess).name, 'data');
        
        P_sig = []; S_sig = [];
        P_MI = [];  S_MI = [];
        P_xcor = [];
        for iC = length(data.SI):-1:1
            P_sig(iC) = data.SI(iC).spatial.place.MI_pval;
            P_MI(iC) = data.SI(iC).spatial.place.MI;
            %                fprintf('Cell %0.0f: MI = %0.3f   pval = %0.3f\n', iC, P_MI(iC), P_sig(iC))
            
            S_sig(iC) = data.SI(iC).spatial.speed.MI_pval;
            S_MI(iC) = data.SI(iC).spatial.speed.MI;
            
            P_xcor(iC) = data.SI(iC).spatial.place.split.Stability_corr;
        end
        all_P_sigs = [all_P_sigs, P_sig];
        all_P_MI = [all_P_MI, P_MI];
        all_P_xcor = [all_P_xcor, P_xcor];
        
        all_S_sigs = [all_S_sigs, S_sig];
        all_S_MI = [all_S_MI, S_MI];
        
        fprintf('<strong>%s  |  Place Sig =  %.2f%%  |  Speed Sig =  %.2f%% </strong>\n', [data.f_info.fname ' ' data.f_info.date ' ' data.f_info.time ' ' data.f_info.task]  , (sum(P_sig < 0.05)/length(P_sig))*100, (sum(S_sig < 0.05)/length(S_sig))*100)
        
        
        all_sess_P_sigs(iSess) = (sum(P_sig < 0.05)/length(P_sig))*100;
        all_sess_P_sigs_split(iSess) = (sum(((P_sig < 0.05) & (P_xcor >= .5)))/length(P_sig))*100;
        
        all_sess_S_sigs(iSess) = (sum(S_sig < 0.05)/length(S_sig))*100;
        
        if data.behav.width <35 && data.behav.height < 35
            fprintf('Small OF  Width: %.0fcm  x  Height: %.0fcm\n', data.behav.width, data.behav.height)
            all_fname{iSess} = [data.f_info.date ' S OF'];
        else
            fprintf('Large OF  Width: %.0fcm  x  Height: %.0fcm\n', data.behav.width, data.behav.height)
            all_fname{iSess} = [data.f_info.date ' L OF'];
        end
        
        
    end % end sessions
    % collect everything
    All_Data.(sub_list{iSub}).all_P_MI = all_P_MI;
    All_Data.(sub_list{iSub}).all_P_sigs = all_P_sigs;
    All_Data.(sub_list{iSub}).all_P_xcor = all_P_xcor;
    
    All_Data.(sub_list{iSub}).all_S_MI = all_S_MI;
    All_Data.(sub_list{iSub}).all_S_sigs = all_S_sigs;
    
    All_Data.(sub_list{iSub}).all_sess_P_sigs = all_sess_P_sigs;
    All_Data.(sub_list{iSub}).all_sess_P_sigs_split = all_sess_P_sigs_split;
    All_Data.(sub_list{iSub}).all_sess_S_sigs = all_sess_S_sigs;
    
    All_Data.(sub_list{iSub}).all_fname = all_fname;
    
    
end % end subjects

%%  go across subjects and plot the % of sig cells per session
sub_list = fieldnames(All_Data); % get the list of subjects

% generate plots across subjects

for iSub = 1:length(sub_list)
    figure(iSub)
    title(sub_list{iSub})
    subplot(6,4,[1:3 5:7 9:11 13:15 17:19])
    c_ord = linspecer(2);
    bh = bar([All_Data.(sub_list{iSub}).all_sess_P_sigs' ; All_Data.(sub_list{iSub}).all_sess_S_sigs']');
    
    bh(1).FaceColor = c_ord(1,:);
    bh(2).FaceColor = c_ord(2,:);
    
    ylabel('% of sig modulated cells')
    set(gca,'xtick', 1:length(All_Data.(sub_list{iSub}).all_fname), 'XTickLabel', All_Data.(sub_list{iSub}).all_fname, 'fontsize', 12)
    xtickangle(45)
    y_lim = ylim;
    
    legend({'Place', 'Speed'})
    
    subplot(6,4,[4 8 12 16 20])
    
    % [~,eh] = errorbar_groups([mean(all_P_sigs) ; mean(all_S_sigs)]', err_bars','bar_colors', c_ord, 'bar_width',0.75,'errorbar_width',0.5, 'bar_names',{'Place', 'Speed'});
    
    err_bars = [std(All_Data.(sub_list{iSub}).all_sess_P_sigs) / sqrt( length(All_Data.(sub_list{iSub}).all_sess_P_sigs)) ; std(All_Data.(sub_list{iSub}).all_sess_S_sigs) / sqrt( length(All_Data.(sub_list{iSub}).all_sess_S_sigs))];
    
    yyaxis right
    bh2 = bar(1 ,[mean(All_Data.(sub_list{iSub}).all_sess_P_sigs) ; mean(All_Data.(sub_list{iSub}).all_sess_S_sigs)]');
    bh2(1).FaceColor = c_ord(1,:);
    bh2(2).FaceColor = c_ord(2,:);
    
    hold on
    eh  = errorbar([.85 1.15],[mean(All_Data.(sub_list{iSub}).all_sess_P_sigs) ; mean(All_Data.(sub_list{iSub}).all_sess_S_sigs)]', err_bars','color', 'k', 'LineStyle',  'none');
    % eh.LineStyle
    
    ylim(y_lim)
    ylabel('% of sig modulated cells')
    % set(gca, 'XTickLabel', all_fname)
    set(gca, 'YColor', 'k','xticklabel', 'Mean')
    
    legend({'Place', 'Speed'})
    yyaxis left
    set(gca, 'ytick', [])
    
end

%% Scatter MIs and MI p vals
% remove zero pvals
c_ord = linspecer(2);


for iSub = 1:length(sub_list)
    
    zero_idx = All_Data.(sub_list{iSub}).all_P_sigs == 0;
    sig_idx = All_Data.(sub_list{iSub}).all_P_sigs < 0.05  & ~zero_idx;
    xcor_idx = All_Data.(sub_list{iSub}).all_P_xcor >= .5;
    
    
    figure(1000+iSub)
    title(strrep(sub_list{iSub}, '_', ' '))
    subplot(1,2,1)
    hold on
    scatter(All_Data.(sub_list{iSub}).all_P_MI(~zero_idx), All_Data.(sub_list{iSub}).all_P_sigs(~zero_idx),36, [.8 .8 .8]);
    scatter(All_Data.(sub_list{iSub}).all_P_MI(sig_idx & xcor_idx), All_Data.(sub_list{iSub}).all_P_sigs(sig_idx & xcor_idx), 36, c_ord(1,:), 'filled');
    
    yline(0.05)
    
    xlabel('MI')
    ylabel('MI pval')
    
    subplot(1,2,2)
    hold on
    scatter(All_Data.(sub_list{iSub}).all_P_MI(sig_idx), All_Data.(sub_list{iSub}).all_P_sigs(sig_idx), 36, [.8 .8 .8])
    scatter(All_Data.(sub_list{iSub}).all_P_MI(sig_idx & xcor_idx), All_Data.(sub_list{iSub}).all_P_sigs(sig_idx & xcor_idx), 36, c_ord(1,:), 'filled')
    
    xlabel('MI')
    ylabel('MI pval')
    
end


%% MI and xcor over sessions Scatter


for iSub = 1:length(sub_list) % comment length out and just do iSub = X session.
    cd(PARAMS.inter_dir)
    
    cd(sub_list{iSub}) % go to the subject folder
    
    sess_list = dir('*.mat');
    all_P_sigs = []; all_P_MI = []; all_sess_P_sigs = NaN(length(sess_list),1); all_sess_P_sigs_split = NaN(length(sess_list),1);  all_P_xcor = [];
    all_S_sigs = []; all_S_MI = []; all_sess_S_sigs = NaN(length(sess_list),1);
    all_fname = cell(length(sess_list),1);
    
    for iSess =  1: length(sess_list)
        
        load(sess_list(iSess).name, 'data');
        
        P_sig = []; S_sig = [];
        P_MI = [];  S_MI = [];
        P_xcor = [];
        for iC = length(data.SI):-1:1
            P_sig(iC) = data.SI(iC).spatial.place.MI_pval;
            P_MI(iC) = data.SI(iC).spatial.place.MI;
            %                fprintf('Cell %0.0f: MI = %0.3f   pval = %0.3f\n', iC, P_MI(iC), P_sig(iC))
            
            S_sig(iC) = data.SI(iC).spatial.speed.MI_pval;
            S_MI(iC) = data.SI(iC).spatial.speed.MI;
            
            P_xcor(iC) = data.SI(iC).spatial.place.split.Stability_corr;
        end
        
        Sess_P_MI{iSess, 1} = P_MI(P_sig < 0.05);
        Sess_P_MI{iSess, 2} = P_MI(P_sig >= 0.05);
        
        Sess_P_xcorr{iSess, 1} = P_xcor(P_sig < 0.05);
        Sess_P_xcorr{iSess, 2} = P_xcor(P_sig >= 0.05);
        
        Sess_date{iSess}  =data.f_info.date;
        fprintf('<strong> %s </strong>\n', [data.f_info.fname ' ' data.f_info.date ' ' data.f_info.time ' ' data.f_info.task])
        
    end
    
    % plot using raincloud
    figure(2000+iSub)
    title(strrep(sub_list{iSub}, '_', ' '))
    hs = rm_raincloud(Sess_P_MI, c_ord);
    
    for ii = 1:length(hs.m)
        hs.s{ii,1}.SizeData = 50;
        hs.s{ii,2}.SizeData = 50;
        
        hs.m(ii,1).SizeData = 50;
        hs.m(ii,2).SizeData = 50;
    end
    for ii = 1:length(hs.l)
        hs.l(ii,1).MarkerSize = 3;
        hs.l(ii,2).MarkerSize = 3;
    end
    
    % adjust.  For some reason rm_raincloud plot flips the x and y
    % axes.
    y_ticks = get(gca, 'ytick');
    ylim([y_ticks(1) y_ticks(end)])
    xlabel('MI')
    ylabel('session number')
    legend({'','sig place','','non-sig'}, 'location', 'North')
    set(gca, 'yTickLabel', Sess_date, 'fontsize', 12)
    ytickangle(45)
    
    
    figure(3000+iSub)
    title(strrep(sub_list{iSub}, '_', ' '))
    hs = rm_raincloud(Sess_P_xcorr, c_ord);
    
    for ii = 1:length(hs.m)
        hs.s{ii,1}.SizeData = 50;
        hs.s{ii,2}.SizeData = 50;
        
        hs.m(ii,1).SizeData = 50;
        hs.m(ii,2).SizeData = 50;
    end
    
    for ii = 1:length(hs.l)
        hs.l(ii,1).MarkerSize = 3;
        hs.l(ii,2).MarkerSize = 3;
    end
    
    % adjust.  For some reason rm_raincloud plot flips the x and y
    % axes.
    y_ticks = get(gca, 'ytick');
    ylim([y_ticks(1) y_ticks(end)])
    xlabel('split xcorr')
    ylabel('session number')
    legend({'','sig place','','non-sig'}, 'location', 'North')
    set(gca, 'yTickLabel', Sess_date, 'fontsize', 12)
    ytickangle(45)
    clear Sess_P_MI Sess_P_xcorr Sess_date
    
end

%%




%% ignore below


%%  simple version
figure(1010)
subplot(6,4,[1:3 5:7 9:11 13:15 17:19])
c_ord = linspecer(2);
bh = bar([all_sess_P_sigs ; all_sess_S_sigs]');

bh(1).FaceColor = c_ord(1,:);
bh(2).FaceColor = c_ord(2,:);

ylabel('% of sig modulated cells')
set(gca,'xtick', 1:length(all_fname), 'XTickLabel', all_fname, 'fontsize', 12)
xtickangle(45)
y_lim = ylim;

legend({'Place', 'Speed'})

subplot(6,4,[4 8 12 16 20])

% [~,eh] = errorbar_groups([mean(all_P_sigs) ; mean(all_S_sigs)]', err_bars','bar_colors', c_ord, 'bar_width',0.75,'errorbar_width',0.5, 'bar_names',{'Place', 'Speed'});

err_bars = [std(all_sess_P_sigs) / sqrt( length(all_sess_P_sigs)) ; std(all_sess_S_sigs) / sqrt( length(all_sess_S_sigs))];

yyaxis right
bh2 = bar(1 ,[mean(all_sess_P_sigs) ; mean(all_sess_S_sigs)]');
bh2(1).FaceColor = c_ord(1,:);
bh2(2).FaceColor = c_ord(2,:);

hold on
eh  = errorbar([.85 1.15],[mean(all_sess_P_sigs) ; mean(all_sess_S_sigs)]', err_bars','color', 'k', 'LineStyle',  'none');
% eh.LineStyle

ylim(y_lim)
ylabel('% of sig modulated cells')
% set(gca, 'XTickLabel', all_fname)
set(gca, 'YColor', 'k','xticklabel', 'Mean')

legend({'Place', 'Speed'})
yyaxis left
set(gca, 'ytick', [])


%% Scatter MIs and MI p vals
% remove zero pvals
c_ord = linspecer(2);

zero_idx = all_P_sigs == 0;
sig_idx = all_P_sigs < 0.05  & ~zero_idx;
xcor_idx = all_P_xcor >= .5;


figure(1000)
subplot(1,2,1)
hold on
scatter(all_P_MI(~zero_idx), all_P_sigs(~zero_idx),36, [.8 .8 .8]);
scatter(all_P_MI(sig_idx & xcor_idx), all_P_sigs(sig_idx & xcor_idx), 36, c_ord(1,:), 'filled');

yline(0.05)

xlabel('MI')
ylabel('MI pval')

subplot(1,2,2)
hold on
scatter(all_P_MI(sig_idx), all_P_sigs(sig_idx), 36, [.8 .8 .8])
scatter(all_P_MI(sig_idx & xcor_idx), all_P_sigs(sig_idx & xcor_idx), 36, c_ord(1,:), 'filled')

xlabel('MI')
ylabel('MI pval')



%%

%             if (iSess == 2 && iTask ==2) || (iSess == 1 && iTask ==1)
%                 continue
%             end
%
%             %% run the screening script
%             These_cells{iSess, iTask} = Spatial_screener_info(cfg, f_info);
%
%             %% check the sig for spatial metrics
%             %         close all
%             fill_space = repmat(' ',1, 30 - length(f_info.fname));
%             fprintf(['<strong>%s</strong>:' fill_space '\n'], f_info.fname)
%             [~, sig_cells] = MS_get_sig_cells(These_cells{iSess, iTask}, 0.01);
%             %
%
%             % spatial/place plots
%             place_sig = find(sig_cells(:,1));
%             if ~isempty(place_sig)
%                 if strcmp(f_info.task, 'LT')
%                     MS_plot_spatial_cell_1D(These_cells{iSess, iTask},place_sig')
%                 else
%                     MS_plot_spatial_cell(These_cells{iSess, iTask},place_sig')
%                 end
%             end
%
%
%             % speed plots
%             speed_sig = find(sig_cells(:,2));
%             if ~isempty(speed_sig)
%                 MS_plot_movement_cell_1D(These_cells{iSess, iTask},speed_sig', 'speed')
%             end
%
%
%             % accel plots
%             accel_sig = find(sig_cells(:,3));
%             if ~isempty(accel_sig)
%                 MS_plot_movement_cell_1D(These_cells{iSess, iTask},speed_sig', 'accel')
%             end
%
%             % cell summary (everything)
%             %         MS_plot_cell(These_cells{iSess, iTask},
%
%


%         end % end tasks
%
%     end %sessions
%
% end% subjects

