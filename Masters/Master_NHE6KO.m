%% MASTER NHE6KO Analyses
%
%   This script will take in long LFP recordings from HC, PFC/S1, and EMG
%   data and perform initial preprocessing (loading, filtering), Mannual
%   sleep state screening, and basic statistics (percentage in
%   wake/SWS/REM, group comparisons of wake/sleep states across time,
%   mean wake/sleep episode durations and number of transitions).
%       These analyses are based on this: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130177
%
%
%   Requirements:
%       - CEH2 codebase: https://github.com/ecarmichael/CEH2
%       -
%%   Initialize paths and set up system parameters
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
    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/'; % raw data location.
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
else
    PARAMS.raw_data_dir = 'J:\Williams_Lab\NHE6KO\Raw_data'; % where to find the raw data
    PARAMS.inter_dir = 'J:\Williams_Lab\NHE6KO\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\NHE6KO\Stats'; % where to put the statistical output .txt
    PARAMS.fig_dir = 'J:\Williams_Lab\NHE6KO\Figs'; % where to store figures.  Best to have subfolders and use dates.
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

% subject specific parameters.  Used for loading specific channels on a subject by subject basis.
PARAMS.Subjects.M02.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M02.LFP_Chan = 'CSC50.ncs';  % best LFP channel.
PARAMS.Subjects.M02.EMG_Chan = 'CSC63.ncs';  % best EMG channel.
PARAMS.Subjects.M02.genotype = 'wt';

PARAMS.Subjects.M05.dir_name = 'NHE6KO_05'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M05.LFP_Chan = 'CSC1.ncs';  % best LFP channel.
PARAMS.Subjects.M05.EMG_Chan = 'CSC15.ncs';  % best EMG channel.
PARAMS.Subjects.M05.genotype = 'wt';

% K/O
% PARAMS.Subjects.M11.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
% PARAMS.Subjects.M11.LFP_Chan = 'CSC39.ncs';  % best LFP channel.
% PARAMS.Subjects.M11.EMG_Chan = 'CSC47.ncs';  % best EMG channel.
% PARAMS.Subjects.M11.genotype = 'ko';

PARAMS.Subjects.M12.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M12.LFP_Chan = 'CSC7.ncs';  % best LFP channel.
PARAMS.Subjects.M12.EMG_Chan = 'CSC16.ncs';  % best EMG channel.
PARAMS.Subjects.M12.genotype = 'ko';

PARAMS.Subjects.M13.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M13.LFP_Chan = 'CSC17.ncs';  % best LFP channel.
PARAMS.Subjects.M13.EMG_Chan = 'CSC32.ncs';  % best EMG channel.
PARAMS.Subjects.M13.genotype = 'ko';

Subjects = fieldnames(PARAMS.Subjects);

%% load and append data across multiple recording blocks.  To prevent buffer errors and corrupted data the 72hr recording period was broken down into several blocks of 8-11 hours (30hrs for mouse #05).


% NH_preprocess(PARAMS.raw_data_dir, PARAMS.inter_dir); % run this script once to load, preprocess, and save all the 'good' channels for each subject to the PARAMS.inter_dir.


%% Classify the sleep states by visual inspection.

for iSub = 1:length(Subjects)
    % load the data for each subject
    load([PARAMS.inter_dir filesep Subjects{iSub} '_raw_data.mat']);
    for iRec = 1:length(this_csc)
        emg_h = abs(hilbert(this_csc{iRec}.data(2,:))); % get the emg power for plotting.
        
        
        cfg_sleep = [];
        cfg_sleep.tvec_range = [0 10];  % number of seconds per window.
        cfg_sleep.emg_range = [min(emg_h) mean(emg_h) + std(emg_h)*5]; % default, should be based on real data.
        cfg_sleep.emg_chan = 1; % emg channel.  Can be empty.
        cfg_sleep.lfp_chans = 1; % lfp channels to be plotted can be empty. Can be 1 or more, but best to keep it less than 3. should be rows in csc.data.
        cfg_sleep.state_name = {'Wake',       'SWS',       'REM',    'Quiescence','Transition','micro',  'Redo',     'Exit'}; %
        cfg_sleep.state_keys = {'rightarrow','uparrow', 'downarrow', 'leftarrow', 'numpad0',   'numpad1' 'backspace','backquote' }; % which key to press for each state_val
        
        % work in hour long blocks to speed things up and prevent memory
        % issues.
        blocks = 1:(3600*this_csc{iRec}.cfg.hdr{1}.SamplingFrequency):length(this_csc{iRec}.tvec);
        score = cell(1,length(blocks));
        for iB  = 1:length(blocks)
            fprintf('Processing Block %d / %d...\n', iB, length(blocks))
            if iB == 1
                score{iB} = MS_Sleep_score_UI(cfg_sleep, this_csc{iRec}.tvec(1:blocks(iB+1)-1), this_csc{iRec}.data(1,1:blocks(iB+1)-1), emg_h(1:blocks(iB+1)-1));
            elseif iB == length(blocks)
                score{iB} = MS_Sleep_score_UI(cfg_sleep, this_csc{iRec}.tvec(blocks(iB):end), this_csc{iRec}.data(1,blocks(iB):end), emg_h(blocks(iB):end));
            else
                score{iB} = MS_Sleep_score_UI(cfg_sleep, this_csc{iRec}.tvec(blocks(iB):blocks(iB+1)-1), this_csc{iRec}.data(1,blocks(iB):blocks(iB+1)-1), emg_h(blocks(iB):blocks(iB+1)-1));
            end
        end % blocks
        
        %         score{iB} = 1*ones(length(this_csc{iRec}.tvec(blocks(iB):end)),1)
        % catch cases where the sleep score will add in too many points to
        % match the range.
        if length(score{end}) ~= length(this_csc{iRec}.tvec(blocks(iB):end))
            %             score{end} =
        end
        
        
        % save the scores
        sleep_score{iRec}.tvec = this_csc{iRec}.tvec(1:length(cat(1,score{:}))); % catch cases where you don't want all the recording blocks (M05)
        sleep_score{iRec}.score = score;
        save([PARAMS.inter_dir filesep Subjects{iSub} '_sleep_data.mat'], 'sleep_score', '-v7.3')
        
        
        % put the score blocks back together
        this_csc{iRec}.score = cat(1,score{:});
        clear emg_h
    end % recordings
end % subjects

%% temporary fix for OG code recordings

% for iS = 1:length(score)
%     if iS == 1
%         score{iS} = score{iS}(1:end-1);
%     elseif iS == length(score)
%         score{iS} = score{iS}(2:end);
%     else
% score{iS} = score{iS}(2:end-1);
%     end
% end
%% append the data from 1 hour scoring blocks into a continous vector with NaN padding between actual recordings. 

% normalize the time across recording blocks
t_start = this_csc{1}.tvec(1);
Z0 = datetime(datestr('04-Sep-2020 08:00:00', 'dd-mmm-yyyy HH:MM:SS'));
clock_start = seconds(this_csc{1}.tstart - Z0);

all_tvec = [];
all_score = [];
all_data = [];
Fs = mode(diff(sleep_score{1}.tvec));
for iRec = 1:length(sleep_score)
    
    % pad with NaNs between recordings;
    if iRec ~= length(sleep_score)
        N_pad_tvec = sleep_score{iRec}.tvec(end):Fs:sleep_score{iRec+1}.tvec(1);
        N_pad_tvec = N_pad_tvec(2:end); % correct for first and last samples in recording block as close as possible
        N_pad = NaN(1,length(N_pad_tvec));
    else
        N_pad_tvec = [];
        N_pad = [];
    end
    
    % tvec with NaN padding between recording blocks
    padded_tvec = [sleep_score{iRec}.tvec;  N_pad_tvec'];
    % collect the tvecs
    all_tvec = [all_tvec; (padded_tvec - t_start)+ double(clock_start)];
    
    %% get all the data with padding
    all_data = [all_data, [this_csc{iRec}.data(1,1:length(this_csc{iRec}.score)), N_pad]];
    
    % collect the score (should be in cells based on 1 hour blocks.
    if iscell(sleep_score{iRec}.score)
        padded_data = [cat(1,sleep_score{iRec}.score{:}); N_pad'];
        all_score = [all_score; padded_data];
    else % odd case which has it as a vector.
        padded_data = [sleep_score{iRec}.score; N_pad']; % add in the end padding
        all_score = [all_score; padded_data];
    end    
    clear padded_tvec paddec_data
end

save([PARAMS.inter_dir filesep Subjects{iSub} '_all_score.mat'], 'all_score', '-v7.3');
save([PARAMS.inter_dir filesep Subjects{iSub} '_all_data.mat'], 'all_data', '-v7.3');
save([PARAMS.inter_dir filesep Subjects{iSub} '_all_tvec.mat'], 'all_tvec', '-v7.3');


%% get the state block durations and transitions

    [start_idx, end_idx, tran_states] = NH_extract_epochs(all_score);
    Fs = mode(diff(all_tvec));
    for iState = length(start_idx):-1:1
       event_len{iState} = round((end_idx{iState} - start_idx{iState}) * Fs); 
       fprintf('State <strong>%d</strong> mean = %0.2fs median = %0.0fs +/- %2.2fs\n', iState, mean(event_len{iState}),median(event_len{iState}), std(event_len{iState})/sqrt(length(event_len{iState})))
    end

%% break into 24 periods
% day_secs = (ceil(clock_start/3600))*3600:(24*3600):all_tvec(end);
max_hrs = (floor((all_tvec(end) - all_tvec(1))/3600)-1)*3600;
hrs = (ceil(clock_start/3600))*3600:3600:(47+ceil(clock_start/3600))*3600; % correct for offset and number of scored hours.
hrs_actual = hrs;
while sum(hrs_actual>=24*3600) >0 % correct for times >24hrs, 48hrs, ...
    hrs_actual(hrs_actual>=24*3600) = hrs_actual(hrs_actual>=24*3600)-24*3600;
end

% collect mean values for percentage of time in each sleep state.
for iT = length(hrs):-1:1
    this_h_idx = nearest_idx([hrs(iT),hrs(iT)+3600],all_tvec);
    disp(hrs(iT)/3600)
    WAKE(iT) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 1)/length(all_score(this_h_idx(1):this_h_idx(2)));
    NREM(iT) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 2)/length(all_score(this_h_idx(1):this_h_idx(2)));
    REM(iT) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 3)/length(all_score(this_h_idx(1):this_h_idx(2)));
end

hour_scores.hrs = hrs;
hour_scores.WAKE = WAKE;
hour_scores.NREM = NREM;
hour_scores.REM = REM;



% split into n x 24 array
hrs_out = circshift(hrs/3600,ceil(clock_start/3600));
WAKE_out = mean(reshape(circshift(WAKE,ceil(clock_start/3600)),24,2),2)';
NREM_out = mean(reshape(circshift(NREM,ceil(clock_start/3600)),24,2),2)';
REM_out = mean(reshape(circshift(REM,ceil(clock_start/3600)),24,2),2)';

%% break out the sleep states and get the normalized PSD
win = 2048/2;

% WAKE
all_wake = all_data;
all_wake(all_score ~= 1) = [];

[WAKE_ppx, WAKE_F] = pwelch(all_wake, hanning(win), win/2, 2*win, 1/Fs);
Norm_idx = nearest_idx([0.125 100],WAKE_F); % should be the same across WAKE/NREM/REM given the windows but recomputed to be safe.

% NREM
all_nrem = all_data;
all_nrem(all_score ~= 2) = [];

[NREM_ppx, NREM_F] = pwelch(all_nrem, hanning(win), win/2, 2*win, 1/Fs);
Norm_idx = nearest_idx([0.125 100],NREM_F);

% REM
all_rem = all_data;
all_rem(all_score ~= 3) = [];

[REM_ppx, REM_F] = pwelch(all_rem, hanning(win), win/2, 2*win, 1/Fs);
Norm_idx = nearest_idx([0.125 100],REM_F);

hour_scores.PSD.WAKE_F = WAKE_F;
hour_scores.PSD.WAKE_ppx = WAKE_ppx;
hour_scores.PSD.NREM_F = NREM_F;
hour_scores.PSD.NREM_ppx = NREM_ppx;
hour_scores.PSD.REM_F = REM_F;
hour_scores.PSD.REM_ppx = REM_ppx;

save([PARAMS.inter_dir filesep Subjects{iSub} '_hour_data.mat'], 'hour_scores', '-v7.3')


%% make some plots
c_ord = linspecer(9);
% tick labels with Z and real time
XTickString{1} = '$$\begin{array}{c} Z0\\8:00\\ \end{array}$$';
XTickString{2} = '$$\begin{array}{c} Z6\\14:00\\ \end{array}$$';
XTickString{3} = '$$\begin{array}{c} Z12\\20:00\\ \end{array}$$';
XTickString{4} = '$$\begin{array}{c} Z18\\2:00\\ \end{array}$$';

figure(1001)
subplot(4,4,1:3)
hold on
rectangle('position', [1 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(1:24, WAKE_out*100,'-o', 'color', c_ord(1,:))
xlabel('time from light onset')
ylabel('% wake');
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
ylim([0 100]); xlim([1 24]);

% bar plot for light vs dark
subplot(4,4,4)
groups = [ones(1,12),2*ones(1,12)];
boxplot(WAKE_out*100, groups);
set(gca, 'xticklabel', {'light', 'dark'})
h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),c_ord(5,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(1,:),'FaceAlpha',.5);


subplot(4,4,5:7)
hold on
rectangle('position', [1 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(1:24, NREM_out*100,'-o', 'color', c_ord(3,:))
xlabel('time from light onset')
ylabel('% NREM');
set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex');
ylim([0 100]); xlim([1 24]);

% bar plot for light vs dark
subplot(4,4,8)
groups = [ones(1,12),2*ones(1,12)];
boxplot(NREM_out*100, groups);
set(gca, 'xticklabel', {'light', 'dark'})
h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),c_ord(5,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(3,:),'FaceAlpha',.5);

subplot(4,4,9:11)
hold on
rectangle('position', [1 0.1 12 20], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(1:24, REM_out*100,'-o', 'color', c_ord(2,:))
xlabel('time from light onset')
ylabel('% REM');
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:2:24, 'xticklabel', {'Z0','Z2', 'Z4', 'Z6','Z8', 'Z10', 'Z12','Z14', 'Z16','Z18','Z20', 'Z22', 'Z24'})
ylim([0 20]); xlim([1 24]);

subplot(4,4,12)
groups = [ones(1,12),2*ones(1,12)];
boxplot(REM_out*100, groups);
set(gca, 'xticklabel', {'light', 'dark'})
h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),c_ord(5,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(2,:),'FaceAlpha',.5);





%% add the PSD plots
figure(1001)
subplot(4,4,13)
plot(WAKE_F, 10*log10(WAKE_ppx), 'color', c_ord(1,:), 'linewidth', 2);
% plot(WAKE_F, ((WAKE_ppx))/sum((WAKE_ppx(Norm_idx))), 'color', c_ord(1,:), 'linewidth', 2);  % using normalization method from Colby-Milley et al. 2015
xlim([0 20])
ylabel('WAKE power')
xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
% subplot(4,4,13)
% axes('Position',[.7 .6 .2 .3])
% box on
% plot(WAKE_F, ((WAKE_ppx))/sum((WAKE_ppx(Norm_idx))), 'color', c_ord(3,:), 'linewidth', 3);
% xlim([20 50])

subplot(4,4,14)
plot(NREM_F, 10*log10(NREM_ppx), 'color', c_ord(3,:), 'linewidth', 2);
% plot(NREM_F, ((NREM_ppx))/sum((NREM_ppx(Norm_idx))), 'color', c_ord(3,:), 'linewidth', 2);
xlim([0 20])
ylabel('NREM power')
xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
% axes('Position',[.7 .6 .2 .3])
% box on
% plot(NREM_F, ((NREM_ppx))/sum((NREM_ppx(Norm_idx))), 'color', c_ord(3,:), 'linewidth', 3);
% xlim([20 50])


subplot(4,4,15)
plot(REM_F, 10*log10(REM_ppx), 'color', c_ord(2,:), 'linewidth', 2);
% plot(REM_F, ((REM_ppx))/sum((REM_ppx(Norm_idx))), 'color', c_ord(2,:), 'linewidth', 2);
xlim([0 20])
ylabel('REM power')
xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
% axes('Position',[.7 .6 .2 .3])
% box on
% plot(REM_F, ((REM_ppx))/sum((REM_ppx(Norm_idx))), 'color', c_ord(3,:), 'linewidth', 3);
% xlim([20 50])


% add in a spare plot with subject info

subplot(4,4,16)
text(0, .75, ['Subject: ' Subjects{iSub}])
text(0, .5, ['Genotype: ' PARAMS.Subjects.(Subjects{iSub}).genotype])
text(0, 0.25, [num2str(length(hour_scores.hrs)/24) 'x24hrs'])
axis off

pos = get(gcf, 'position');
set(gcf, 'position', [pos(1) pos(2)*.6 pos(3)*1.4 pos(4)*1.4]);

saveas(gcf, [PARAMS.inter_dir filesep Subjects{iSub} 'summary.png']);
saveas(gcf, [PARAMS.inter_dir filesep Subjects{iSub} 'summary.fig']);

%%  Cross subject data table
% types = {'ko', 'wt
Sub_id = [];
Gene_id = {};
Z_val = [];  % time in Zeitgeber time Z0 = lights turn on. 
LD_val = {}; % light or dak cycle value for each hour. 
wake_val = [];
NREM_val = [];
REM_val = [];

for iSub = 1:length(Subjects)
%     if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}) = load([PARAMS.inter_dir filesep Subjects{iSub} '_hour_data.mat']);
%     elseif contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
%         Sub_scores.wt.(Subjects{iSub}) = load([PARAMS.inter_dir filesep Subjects{iSub} '_hour_data.mat']);
%     end
%     
    wake_val = [wake_val, circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.WAKE,ceil(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.hrs(1)/3600))];
    NREM_val = [NREM_val,circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.NREM,ceil(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.hrs(1)/3600))];
    REM_val = [REM_val,circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.REM,ceil(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.hrs(1)/3600))];
    
    Sub_id = [Sub_id, repmat(str2double(Subjects{iSub}(2:end)),1,length(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.hrs))];
    gene_temp = cell(1,length(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.hrs)); 
    gene_temp(:) = {PARAMS.Subjects.(Subjects{iSub}).genotype};
    Gene_id = [Gene_id, gene_temp]; 
    Z_val = [Z_val, [0:23 0:23]]; 
        LD_temp = cell(1,length(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.hrs)); 
    LD_temp(1:12) = {'Light'}; LD_temp(13:24) = {'Dark'}; LD_temp(25:36) = {'Light'}; LD_temp(37:48) = {'Dark'}; 
%     LD_binary = [ones(1,12) zeros(1,12) ones(1,12) zeros(1,12)]; 
    LD_val = [LD_val, LD_temp]; 
end


NHE6_tbl = table(nominal(Sub_id), nominal(Gene_id), ordinal(Z_val), nominal(LD_val),wake_val, NREM_val, REM_val, 'VariableNames',{'Subject', 'Genotype', 'Ztime','Light_Dark', 'Wake', 'NREM', 'REM'}); 

%%

c_ord = linspecer(9);
% tick labels with Z and real time
XTickString{1} = '$$\begin{array}{c} Z0\\8:00\\ \end{array}$$';
XTickString{2} = '$$\begin{array}{c} Z6\\14:00\\ \end{array}$$';
XTickString{3} = '$$\begin{array}{c} Z12\\20:00\\ \end{array}$$';
XTickString{4} = '$$\begin{array}{c} Z18\\2:00\\ \end{array}$$';

figure(1001)
subplot(5,4,1:3)
hold on
rectangle('position', [1 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [13 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0]);

ko_mean = []; wt_mean = [];
for iSub = 1:length(Subjects)
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        this_val =Sub_scores.ko.(Subjects{iSub}).hour_scores.WAKE;
        this_time = Sub_scores.ko.(Subjects{iSub}).hour_scores.hrs/3600;
        p(iSub) = plot(1:24, mean(reshape(circshift(this_val,ceil(this_time(1))),24,2),2)'*100,'-.', 'color',[c_ord(1,:) .2], 'MarkerEdgeColor',c_ord(1,:));
        ko_mean = [ko_mean; this_val];
    end
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        this_val =Sub_scores.wt.(Subjects{iSub}).hour_scores.WAKE;
        this_time = Sub_scores.wt.(Subjects{iSub}).hour_scores.hrs/3600;
        p(iSub) = plot(1:24, mean(reshape(circshift(this_val,ceil(this_time(1))),24,2),2)'*100,'--', 'color',[0 0 0 .2], 'MarkerEdgeColor',[0 0 0]);
        wt_mean = [wt_mean; this_val];
    end
    
    
end
p(iSub+1) = plot(1:24, median(reshape(circshift(mean(wt_mean),ceil(this_time(1))),24,2),2)'*100,'-*', 'color','k', 'MarkerEdgeColor','k', 'linewidth',1.5 );
p(iSub+2) = plot(1:24, median(reshape(circshift(mean(ko_mean),ceil(this_time(1))),24,2),2)'*100,'-x', 'color',c_ord(1,:), 'MarkerEdgeColor',c_ord(1,:), 'linewidth', 1.5 );
% plot(1:24, WAKE_out*100,'-o', 'color', c_ord(1,:))
% xlabel('time from light onset')
ylabel('% wake');
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
ylim([0 100]); xlim([1 24]);
legend([p(iSub+1) p(iSub+2)],{'wt', 'ko'}, 'box', 'off', 'orientation', 'horizontal');


% bar plot for light vs dark
subplot(5,4,4)
hold on
rectangle('position', [0, 0, 3, 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [3, 0, 3, 100], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0]);

% data = [NHE6_tbl.Wake(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Light')*100; NHE6_tbl.Wake(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Light')*100;...
% %     NaN(size(NHE6_tbl.Wake(NHE6_tbl.Genotype == 'ko'& NHE6_tbl.Light_Dark == 'Light'))); ...
%     NHE6_tbl.Wake(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Dark')*100; NHE6_tbl.Wake(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Dark')*100]; 

data = [reshape(circshift(wt_mean, ceil(this_time(1)),2), 1, numel(wt_mean)),reshape(circshift(ko_mean, ceil(this_time(1)),2), 1, numel(ko_mean))] ; 
groups = [repmat([ones(1,12) 3*ones(1,12)],1,2),repmat([ones(1,12) 3*ones(1,12)],1,2), repmat([2*ones(1,12) 4*ones(1,12)],1,2),repmat([2*ones(1,12) 4*ones(1,12)],1,2)]; % groups ([wtL1-12 wtL25-36] [koL1-12 koL25-36] [wtD1-12 wtD25-36] [koD1-12 koD25-36]
% groups = [ones(1,length(data)), 2*ones(1,length(data)),  4*ones(1,length(data)), 5*ones(1,length(data))]; 
boxplot(data'*100, groups, 'Notch','on', 'Positions', [1,2,4,5])
% groups = [ones(1,12),2*ones(1,12)];
% boxplot(WAKE_out*100, groups);
set(gca, 'xticklabel', {'w/t', 'k/o','w/t', 'k/o'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(1,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),c_ord(1,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),'k','FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),'k','FaceAlpha',.5);
ylim([0 100])


% same for NREM
subplot(5,4,5:7)
hold on
rectangle('position', [1 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [13 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0])
ko_mean = []; wt_mean = [];
for iSub = 1:length(Subjects)
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        this_val =Sub_scores.ko.(Subjects{iSub}).hour_scores.NREM;
        this_time = Sub_scores.ko.(Subjects{iSub}).hour_scores.hrs/3600;
        p(iSub) = plot(1:24, mean(reshape(circshift(this_val,ceil(this_time(1))),24,2),2)'*100,'-.', 'color',[c_ord(3,:) .2], 'MarkerEdgeColor',c_ord(3,:));
        ko_mean = [ko_mean; this_val];
    end
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        this_val =Sub_scores.wt.(Subjects{iSub}).hour_scores.NREM;
        this_time = Sub_scores.wt.(Subjects{iSub}).hour_scores.hrs/3600;
        p(iSub) = plot(1:24, mean(reshape(circshift(this_val,ceil(this_time(1))),24,2),2)'*100,'--', 'color',[0 0 0 .2], 'MarkerEdgeColor',[0 0 0]);
        wt_mean = [wt_mean; this_val];
    end
end
p(iSub+1) = plot(1:24, median(reshape(circshift(mean(wt_mean),ceil(this_time(1))),24,2),2)'*100,'-*', 'color','k', 'MarkerEdgeColor','k', 'linewidth',1.5 );
p(iSub+2) = plot(1:24, median(reshape(circshift(mean(ko_mean),ceil(this_time(1))),24,2),2)'*100,'-x', 'color',c_ord(3,:), 'MarkerEdgeColor',c_ord(3,:), 'linewidth', 1.5 );
% plot(1:24, NREM_out*100,'-o', 'color', c_ord(3,:))
% xlabel('time from light onset')
ylabel('% NREM');
set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex');
ylim([0 100]); xlim([1 24]);
legend([p(iSub+1) p(iSub+2)],{'wt', 'ko'}, 'box', 'off', 'orientation', 'horizontal');

subplot(5,4,8)
hold on
rectangle('position', [0, 0, 3, 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [3, 0, 3, 100], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0]);

% data = [NHE6_tbl.NREM(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Light')*100; NHE6_tbl.NREM(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Light')*100;...
% %     NaN(size(NHE6_tbl.Wake(NHE6_tbl.Genotype == 'ko'& NHE6_tbl.Light_Dark == 'Light'))); ...
%     NHE6_tbl.NREM(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Dark')*100; NHE6_tbl.NREM(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Dark')*100]; 
% groups = [ones(1,length(data)), 2*ones(1,length(data)),  4*ones(1,length(data)), 5*ones(1,length(data))]; 
data = [reshape(circshift(wt_mean, ceil(this_time(1)),2), 1, numel(wt_mean)),reshape(circshift(ko_mean, ceil(this_time(1)),2), 1, numel(ko_mean))] ; 
groups = [repmat([ones(1,12) 3*ones(1,12)],1,2),repmat([ones(1,12) 3*ones(1,12)],1,2), repmat([2*ones(1,12) 4*ones(1,12)],1,2),repmat([2*ones(1,12) 4*ones(1,12)],1,2)]; % groups ([wtL1-12 wtL25-36] [koL1-12 koL25-36] [wtD1-12 wtD25-36] [koD1-12 koD25-36]
boxplot(data'*100, groups, 'Notch','on', 'Positions', [1,2,4,5])
% groups = [ones(1,12),2*ones(1,12)];
% boxplot(WAKE_out*100, groups);
set(gca, 'xticklabel', {'w/t', 'k/o','w/t', 'k/o'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(3,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),c_ord(3,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),'k','FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),'k','FaceAlpha',.5);
ylim([0 100])


subplot(5,4,9:11)
hold on
rectangle('position', [1 0.1 12 20], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [13 0.1 12 20], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0])
ko_mean = []; wt_mean = [];
for iSub = 1:length(Subjects)
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        this_val =Sub_scores.ko.(Subjects{iSub}).hour_scores.REM;
        this_time = Sub_scores.ko.(Subjects{iSub}).hour_scores.hrs/3600;
        p(iSub) = plot(1:24, mean(reshape(circshift(this_val,ceil(this_time(1))),24,2),2)'*100,'-.', 'color',[c_ord(2,:) .2], 'MarkerEdgeColor',c_ord(2,:));
        ko_mean = [ko_mean; this_val];
    end
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        this_val =Sub_scores.wt.(Subjects{iSub}).hour_scores.REM;
        this_time = Sub_scores.wt.(Subjects{iSub}).hour_scores.hrs/3600;
        p(iSub) = plot(1:24, mean(reshape(circshift(this_val,ceil(this_time(1))),24,2),2)'*100,'--', 'color',[0 0 0 .2], 'MarkerEdgeColor',[0 0 0]);
        wt_mean = [wt_mean; this_val];
    end
end
p(iSub+1) = plot(1:24, median(reshape(circshift(mean(wt_mean),ceil(this_time(1))),24,2),2)'*100,'-*', 'color','k', 'MarkerEdgeColor','k', 'linewidth',1.5 );
p(iSub+2) = plot(1:24, median(reshape(circshift(mean(ko_mean),ceil(this_time(1))),24,2),2)'*100,'-x', 'color',c_ord(2,:), 'MarkerEdgeColor',c_ord(2,:), 'linewidth', 1.5 );
% plot(1:24, REM_out*100,'-o', 'color', c_ord(2,:))
% xlabel('time from light onset')
ylabel('% REM');
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:2:24, 'xticklabel', {'Z0','Z2', 'Z4', 'Z6','Z8', 'Z10', 'Z12','Z14', 'Z16','Z18','Z20', 'Z22', 'Z24'})
ylim([0 20]); xlim([1 24]);
legend([p(iSub+1) p(iSub+2)],{'wt', 'ko'}, 'box', 'off', 'orientation', 'horizontal');

subplot(5,4,12)
hold on
rectangle('position', [0, 0, 3, 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
rectangle('position', [3, 0, 3, 100], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0])

% data = [NHE6_tbl.REM(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Light')*100; NHE6_tbl.REM(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Light')*100;...
% %     NaN(size(NHE6_tbl.Wake(NHE6_tbl.Genotype == 'ko'& NHE6_tbl.Light_Dark == 'Light'))); ...
%     NHE6_tbl.REM(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Dark')*100; NHE6_tbl.REM(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Dark')*100]; 
% % groups = [ones(1,length(data)), 2*ones(1,length(data)),  4*ones(1,length(data)), 5*ones(1,length(data))]; 
% boxplot(data', 'Notch','on', 'Positions', [1,2,4,5])
data = [reshape(circshift(wt_mean, ceil(this_time(1)),2), 1, numel(wt_mean)),reshape(circshift(ko_mean, ceil(this_time(1)),2), 1, numel(ko_mean))] ; 
groups = [repmat([ones(1,12) 3*ones(1,12)],1,2),repmat([ones(1,12) 3*ones(1,12)],1,2), repmat([2*ones(1,12) 4*ones(1,12)],1,2),repmat([2*ones(1,12) 4*ones(1,12)],1,2)]; % groups ([wtL1-12 wtL25-36] [koL1-12 koL25-36] [wtD1-12 wtD25-36] [koD1-12 koD25-36]
boxplot(data'*100, groups, 'Notch','on', 'Positions', [1,2,4,5])
% groups = [ones(1,12),2*ones(1,12)];
% boxplot(WAKE_out*100, groups);
set(gca, 'xticklabel', {'w/t', 'k/o','w/t', 'k/o'})
h = findobj(gca,'Tag','Box');
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(2,:),'FaceAlpha',.5);
patch(get(h(3),'XData'),get(h(3),'YData'),c_ord(2,:),'FaceAlpha',.5);
patch(get(h(2),'XData'),get(h(2),'YData'),'k','FaceAlpha',.5);
patch(get(h(4),'XData'),get(h(4),'YData'),'k','FaceAlpha',.5);
ylim([0 20])

% add in the mean PSDs
wt_mean_wake_psd = [];
wt_mean_NREM_psd = [];
wt_mean_REM_psd = [];

ko_mean_wake_psd = [];
ko_mean_NREM_psd = [];
ko_mean_REM_psd = [];

%%%%%%%% fix this by normalizing in this loop.  also need to get PSDs saved

for iSub = 1:length(Subjects)
    this_WAKE = Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.PSD.WAKE_ppx;
    this_NREM = Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.PSD.NREM_ppx;
    this_REM = Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.PSD.REM_ppx;
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        F_range = nearest_idx([1.25, 100],Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.PSD.WAKE_F); % worksout to be very close even with very subtle difference sin frequencies.
        wt_mean_wake_psd = [wt_mean_wake_psd, this_WAKE./sum(this_WAKE(F_range(1):F_range(2)))];
        wt_mean_NREM_psd = [wt_mean_NREM_psd, this_NREM./sum(this_NREM(F_range(1):F_range(2)))];
        wt_mean_REM_psd = [wt_mean_REM_psd, this_REM./sum(this_REM(F_range(1):F_range(2)))];
    elseif contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        F_range = nearest_idx([1.25, 100],Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).hour_scores.PSD.WAKE_F); % worksout to be very close even with very subtle difference sin frequencies.
        ko_mean_wake_psd = [ko_mean_wake_psd, this_WAKE./sum(this_WAKE(F_range(1):F_range(2)))];
        ko_mean_NREM_psd = [ko_mean_NREM_psd, this_NREM./sum(this_NREM(F_range(1):F_range(2)))];
        ko_mean_REM_psd = [ko_mean_REM_psd, this_REM./sum(this_REM(F_range(1):F_range(2)))];
    end
end

figure(1001)
subplot(5,4,13)
hold on
plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.WAKE_F, median(wt_mean_wake_psd,2), 'color', [0 0 0 .2], 'linewidth', 2);
% plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.WAKE_F, median(wt_mean_wake_psd,2),'--', 'color', c_ord(1,:), 'linewidth', 2);
h1 = shadedErrorBar(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.WAKE_F, median(wt_mean_wake_psd,2), std(wt_mean_wake_psd,0,2)/sqrt(length(wt_mean_wake_psd)));
h1.mainLine.Color = c_ord(1, :);
h1.mainLine.LineWidth =2;
h1.mainLine.LineStyle ='--';
h1.patch.FaceColor = c_ord(1, :);
h1.patch.EdgeColor = c_ord(1, :);
h1.patch.FaceAlpha = .2;
h1.patch.EdgeAlpha = .2;
xlim([0 20])
% ylabel('normalized power')
% xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
y_val = ylim;
text(11, y_val(2)*.9, 'w/t Wake', 'color', [c_ord(1,:) .2], 'fontweight', 'bold')
% legend({'Wake'}, 'box', 'off','FontSize',10)

subplot(5,4,14)
hold on
plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.NREM_F, median(wt_mean_NREM_psd,2), 'color',[0 0 0 .2], 'linewidth', 2);
% plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.NREM_F, median(wt_mean_NREM_psd,2),'--', 'color',c_ord(3,:), 'linewidth', 2);
h1 = shadedErrorBar(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.NREM_F, median(wt_mean_NREM_psd,2), std(wt_mean_NREM_psd,0,2)/sqrt(length(wt_mean_NREM_psd)));
h1.mainLine.Color = c_ord(3, :);
h1.mainLine.LineWidth =2;
h1.mainLine.LineStyle ='--';
h1.patch.FaceColor = c_ord(3, :);
h1.patch.EdgeColor = c_ord(3, :);
h1.patch.FaceAlpha = .2;
h1.patch.EdgeAlpha = .2;
xlim([0 20])
% ylabel('w/t normalized power')
% xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
% legend({'NREM'}, 'box', 'off','FontSize',10)
y_val = ylim;
text(11,y_val(2)*.9, 'w/t NREM', 'color', [c_ord(3,:) .2], 'fontweight', 'bold')


subplot(5,4,15)
hold on
plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.REM_F, median(wt_mean_REM_psd,2), 'color', [0 0 0 .2], 'linewidth', 2);
% plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.REM_F, median(wt_mean_REM_psd,2),'--', 'color', c_ord(2,:), 'linewidth', 2);
h1 = shadedErrorBar(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.REM_F, median(wt_mean_REM_psd,2), std(wt_mean_REM_psd,0,2)/sqrt(length(wt_mean_REM_psd)));
h1.mainLine.Color = c_ord(2, :);
h1.mainLine.LineWidth =2;
h1.mainLine.LineStyle ='--';
h1.patch.FaceColor = c_ord(2, :);
h1.patch.EdgeColor = c_ord(2, :);
h1.patch.FaceAlpha = .2;
h1.patch.EdgeAlpha = .2;
xlim([0 20])
% ylabel('w/t normalized power')
% xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
y_val = ylim;
text(11, y_val(2)*.9, 'w/t REM', 'color', [c_ord(2,:) .2], 'fontweight', 'bold')
% legend({'REM'}, 'box', 'off','FontSize',10)

subplot(5,4,17)
% plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.WAKE_F, median(ko_mean_wake_psd,2), 'color', c_ord(1,:), 'linewidth', 2);
h1 = shadedErrorBar(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.WAKE_F, median(ko_mean_wake_psd,2), std(ko_mean_wake_psd,0,2)/sqrt(length(ko_mean_wake_psd)));
h1.mainLine.Color = c_ord(1, :);
h1.mainLine.LineWidth =2;
h1.patch.FaceColor = c_ord(1, :);
h1.patch.EdgeColor = c_ord(1, :);
h1.patch.FaceAlpha = .2;
h1.patch.EdgeAlpha = .2;
xlim([0 20])
ylabel('normalized power')
xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
y_val = ylim;
text(11, y_val(2)*.9, 'k/o Wake', 'color', c_ord(1,:), 'fontweight', 'bold')
% legend({'Wake'}, 'box', 'off','FontSize',10)

subplot(5,4,18)
% plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.NREM_F, median(ko_mean_NREM_psd,2), 'color',c_ord(3,:), 'linewidth', 2);
h1 = shadedErrorBar(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.NREM_F, median(ko_mean_NREM_psd,2), std(ko_mean_NREM_psd,0,2)/sqrt(length(ko_mean_NREM_psd)));
h1.mainLine.Color = c_ord(3, :);
h1.mainLine.LineWidth =2;
h1.patch.FaceColor = c_ord(3, :);
h1.patch.EdgeColor = c_ord(3, :);
h1.patch.FaceAlpha = .2;
h1.patch.EdgeAlpha = .2;
xlim([0 20])
% ylabel('normalized power')
% xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
% legend({'NREM'}, 'box', 'off','FontSize',10)
y_val = ylim;
text(11, y_val(2)*.9, 'k/o NREM', 'color', c_ord(3,:), 'fontweight', 'bold')


subplot(5,4,19)
hold on
plot(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.REM_F, median(ko_mean_REM_psd,2), 'color', c_ord(2,:), 'linewidth', 2);
% errorb(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.REM_F, median(ko_mean_REM_psd,2), std(ko_mean_REM_psd,0,2)/sqrt(length(ko_mean_REM_psd))); 
h1 = shadedErrorBar(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.PSD.REM_F, median(ko_mean_REM_psd,2), std(ko_mean_REM_psd,0,2)/sqrt(length(ko_mean_REM_psd)));
h1.mainLine.Color = c_ord(2, :);
h1.mainLine.LineWidth =2;
h1.patch.FaceColor = c_ord(2, :);
h1.patch.EdgeColor = c_ord(2, :);
h1.patch.FaceAlpha = .2;
h1.patch.EdgeAlpha = .2;
xlim([0 20])
% ylabel('normalized power')
% xlabel('frequency (Hz)')
set(gca, 'xtick', 0:5:20);
y_val = ylim;
text(11, y_val(2)*.9, 'k/o REM', 'color', c_ord(2,:), 'fontweight', 'bold')
% legend({'REM'}, 'box', 'off','FontSize',10)


subplot(5,4,16)
text(0, .5, 'Median values across subjects')
text(0, .9, 'Genotype: w/t [M02 M05]')
text(0, 0.1, [num2str(length(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.hrs)/24) ' x 24hrs'])
axis off


subplot(5,4,20)
text(0, .9, 'Genotype: k/o [M12 M13]')
text(0, .5, 'Median values across subjects')
text(0, 0.1, [num2str(length(Sub_scores.(PARAMS.Subjects.(Subjects{1}).genotype).(Subjects{1}).hour_scores.hrs)/24) ' x 24hrs'])
pos = get(gcf, 'position');
axis off

cfg_plot.ft_size = 12;
SetFigure(cfg_plot, gcf);
set(gcf, 'position', [pos(1) pos(2)*.4 pos(3)*1.4 pos(4)*1.4]);
saveas(gcf, [PARAMS.inter_dir filesep 'All_subject_summary.png']);
saveas(gcf, [PARAMS.inter_dir filesep 'All_subject_summary.fig']);

%%  Get the epoch stats across subjects


for iSub = 1:length(Subjects)
    fprintf('Loading Subject: %s...\n', Subjects{iSub})
    load([PARAMS.inter_dir filesep Subjects{iSub} '_sleep_data.mat'])
    score_out = []; 
    for iS =1:length(sleep_score)
        score_out = [score_out, sleep_score{iS}.score];
    end
    [start_idx, end_idx, tran_val] = NH_extract_epochs(cat(1,score_out{:})); 
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% same thing but in minutes  and add in event duration, transition, rate values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval = 10; % in minutes. 
% mins = ((ceil(clock_start/3600))*60):interval:((47*60)+(ceil(clock_start/3600)*60)); % correct for offset and number of scored hours.
mins = (ceil(clock_start/3600))*3600:(interval*60):(48+ceil(clock_start/3600))*3600; 
mins = mins(1:end-1); 
mins_actual = mins;
while sum(mins_actual>=48*3600) >0 % correct for times >24hrs, 48hrs, ...
    mins_actual(mins_actual>=48*3600) = mins_actual(mins_actual>=48*3600)-(48*3600);
end

WAKE = []; NREM = []; REM = []; 
% collect mean values for percentage of time in each sleep state.
for iT = length(mins):-1:1
    this_h_idx = nearest_idx([mins(iT),mins(iT)+interval*60],all_tvec);
    fprintf('Block: %d  df: %.0fsec ...\n', mins(iT)/60, all_tvec(this_h_idx(2)) - all_tvec(this_h_idx(1)))
    WAKE(iT) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 1)/length(all_score(this_h_idx(1):this_h_idx(2)));
    NREM(iT) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 2)/length(all_score(this_h_idx(1):this_h_idx(2)));
    REM(iT) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 3)/length(all_score(this_h_idx(1):this_h_idx(2)));
end



min_scores.mins = mins_actual;
min_scores.WAKE = WAKE;
min_scores.NREM = NREM;
min_scores.REM = REM;
min_scores.clock_start = clock_start; 


% split into n x 24 array
mins_out = circshift(min_scores.mins/60,((ceil(min_scores.clock_start/3600))*(60/interval)));
WAKE_out = circshift(min_scores.WAKE,((ceil(min_scores.clock_start/3600))*(60/interval)));
NREM_out = circshift(min_scores.NREM,((ceil(min_scores.clock_start/3600))*(60/interval)));
REM_out = circshift(min_scores.REM,((ceil(min_scores.clock_start/3600))*(60/interval)));

save([PARAMS.inter_dir filesep Subjects{iSub} '_min_data.mat'], 'min_scores', '-v7.3')

%% make some plots
c_ord = linspecer(9);

groups = [ones(1,length(0:1/(60/interval):12)-1),2*ones(1,length(0:1/(60/interval):12)-1),ones(1,length(0:1/(60/interval):12)-1),2*ones(1,length(0:1/(60/interval):12)-1)];

figure(1002)
subplot(4,4,1:3)
hold on
rectangle('position', [0 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
rectangle('position', [24 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(mins_out/60, WAKE_out*100, 'color', c_ord(1,:))
xlabel('time from light onset')
ylabel('% wake');
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
set(gca, 'xtick', 0:12:48, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
ylim([0 100]); xlim([0 48]);

% bar plot for light vs dark
subplot(4,4,4)
boxplot(WAKE_out*100, groups);
set(gca, 'xticklabel', {'light', 'dark'})
h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),c_ord(5,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(1,:),'FaceAlpha',.5);


subplot(4,4,5:7)
hold on
rectangle('position', [0 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
rectangle('position', [24 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(mins_out/60, NREM_out*100,'color', c_ord(4,:))
xlabel('time from light onset')
ylabel('% NREM');
set(gca, 'xtick', 0:12:48, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex');
ylim([0 100]); xlim([0 48]);

% bar plot for light vs dark
subplot(4,4,8)
boxplot(NREM_out*100, groups);
set(gca, 'xticklabel', {'light', 'dark'})
h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),c_ord(5,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(3,:),'FaceAlpha',.5);

subplot(4,4,9:11)
hold on
rectangle('position', [0 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
rectangle('position', [24 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(mins_out/60, REM_out*100, 'color', c_ord(2,:))
xlabel('time from light onset')
ylabel('% REM');
% set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
set(gca, 'xtick', 1:12:48, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:2:24, 'xticklabel', {'Z0','Z2', 'Z4', 'Z6','Z8', 'Z10', 'Z12','Z14', 'Z16','Z18','Z20', 'Z22', 'Z24'})
ylim([0 40]); xlim([0 48]);

subplot(4,4,12)
boxplot(REM_out*100, groups);
set(gca, 'xticklabel', {'light', 'dark'})
h = findobj(gca,'Tag','Box');
patch(get(h(2),'XData'),get(h(2),'YData'),c_ord(5,:),'FaceAlpha',.5);
patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(2,:),'FaceAlpha',.5);



saveas(gcf, [PARAMS.inter_dir filesep Subjects{iSub} '_' num2str(interval) '_min_summary.png']);
saveas(gcf, [PARAMS.inter_dir filesep Subjects{iSub} '_' num2str(interval) '_min_summary.fig']);


%% 
%%  Cross subject data table using minutes
% types = {'ko', 'wt
Sub_id = [];
Gene_id = {};
Z_val = [];  % time in Zeitgeber time Z0 = lights turn on. 
LD_val = {}; % light or dak cycle value for each hour. 
wake_val = [];
NREM_val = [];
REM_val = [];
min_val = [];
for iSub = 1:length(Subjects)
%     if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}) = load([PARAMS.inter_dir filesep Subjects{iSub} '_min_data.mat']);
%     elseif contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
%         Sub_scores.wt.(Subjects{iSub}) = load([PARAMS.inter_dir filesep Subjects{iSub} '_hour_data.mat']);
%     end
%     
interval = mode(diff(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins/60));  % get the minute interval from the time vector; 
shift_val = ceil((Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.clock_start/3600))*(60/interval); 

    this_min = circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins,shift_val);
    min_val = [min_val, circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins,shift_val)];
    wake_val = [wake_val, circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.WAKE,shift_val)];
    NREM_val = [NREM_val,circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.NREM,shift_val)];
    REM_val = [REM_val,circshift(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.REM,shift_val)];
    
    Sub_id = [Sub_id, repmat(str2double(Subjects{iSub}(2:end)),1,length(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins))];
    gene_temp = cell(1,length(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins)); 
    gene_temp(:) = {PARAMS.Subjects.(Subjects{iSub}).genotype};
    Gene_id = [Gene_id, gene_temp]; 
    Z_val = [Z_val, this_min]; 
        LD_temp = cell(1,length(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins)); 
        L_idx = this_min <= 12*3600 | (this_min >= 24*3600  & this_min <= 36*3600); 
        LD_temp(L_idx) = {'Light'}; 
        LD_temp(~L_idx) = {'Dark'}; 
%     LD_temp(1:interval/12) = {'Light'}; LD_temp(13:24) = {'Dark'}; LD_temp(25:36) = {'Light'}; LD_temp(37:48) = {'Dark'}; 
%     LD_binary = [ones(1,12) zeros(1,12) ones(1,12) zeros(1,12)]; 
    LD_val = [LD_val, LD_temp]; 
end


NHE6_tbl = table(nominal(Sub_id), nominal(Gene_id), ordinal(Z_val), nominal(LD_val),wake_val, NREM_val, REM_val, 'VariableNames',{'Subject', 'Genotype', 'Ztime','Light_Dark', 'Wake', 'NREM', 'REM'}); 


%%
c_ord = linspecer(9);

figure(1001)
subplot(3,4,1:3)
hold on
rectangle('position', [0 0.1 12 100], 'facecolor',[c_ord(5,:), 0.3], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [12 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.1], 'edgecolor', [.02 .02 .02 0]);
rectangle('position', [24 0.1 12 100], 'facecolor',[c_ord(5,:), 0.3], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [36 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.1], 'edgecolor', [.02 .02 .02 0]);
ko_mean = []; wt_mean = [];
for iSub = 1:length(Subjects)
    interval = mode(diff(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins/60));  % get the minute interval from the time vector; 
    shift_val = ceil((Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.clock_start/3600))*(60/interval); 

    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        this_val =circshift(Sub_scores.ko.(Subjects{iSub}).min_scores.WAKE, shift_val);
        this_time = circshift(Sub_scores.ko.(Subjects{iSub}).min_scores.mins/3600, shift_val);
        ko_mean = [ko_mean; this_val];
    end
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        this_val =circshift(Sub_scores.wt.(Subjects{iSub}).min_scores.WAKE, shift_val);
        this_time = circshift(Sub_scores.wt.(Subjects{iSub}).min_scores.mins/3600, shift_val);
        wt_mean = [wt_mean; this_val];
    end
end
% p(iSub+1) = plot(this_time, median(wt_mean)'*100,'-*', 'color','k', 'MarkerEdgeColor','k', 'linewidth',1.5 );
% p(iSub+2) = plot(this_time, median(ko_mean)'*100,'-x', 'color',c_ord(1,:), 'MarkerEdgeColor',c_ord(1,:), 'linewidth', 1.5 );
p(iSub+2) = bar(this_time, median(ko_mean)'*100, 'facecolor', 'k', 'facealpha',0.7);
p(iSub+1) = bar(this_time, median(wt_mean)'*100, 'facecolor', c_ord(1,:),'facealpha', 1);
ylabel('% wake');
set(gca, 'xtick', 0:12:48, 'xticklabel', {'Z0', 'Z12', 'Z24', 'Z36', 'Z48'})
ylim([0 100]); xlim([0 48]);
legend([p(iSub+1) p(iSub+2)],{'wt', 'ko'}, 'box', 'off', 'orientation', 'horizontal');

subplot(3,4,4)
% text(.5, .9, 


% same for NREM
subplot(3,4,5:7)
hold on
rectangle('position', [0 0.1 12 100], 'facecolor',[c_ord(5,:), 0.3], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [12 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.1], 'edgecolor', [.02 .02 .02 0]);
rectangle('position', [24 0.1 12 100], 'facecolor',[c_ord(5,:), 0.3], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [36 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.1], 'edgecolor', [.02 .02 .02 0]);
ko_mean = []; wt_mean = [];
for iSub = 1:length(Subjects)
    interval = mode(diff(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins/60));  % get the minute interval from the time vector; 
    shift_val = ceil((Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.clock_start/3600))*(60/interval); 

    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        this_val =circshift(Sub_scores.ko.(Subjects{iSub}).min_scores.NREM, shift_val);
        this_time = circshift(Sub_scores.ko.(Subjects{iSub}).min_scores.mins/3600, shift_val);
        ko_mean = [ko_mean; this_val];
    end
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        this_val =circshift(Sub_scores.wt.(Subjects{iSub}).min_scores.NREM, shift_val);
        this_time = circshift(Sub_scores.wt.(Subjects{iSub}).min_scores.mins/3600, shift_val);
        wt_mean = [wt_mean; this_val];
    end
end
% p(iSub+1) = plot(this_time, median(wt_mean)'*100,'-*', 'color','k', 'MarkerEdgeColor','k', 'linewidth',1.5 );
% p(iSub+2) = plot(this_time, median(ko_mean)'*100,'-x', 'color',c_ord(3,:), 'MarkerEdgeColor',c_ord(3,:), 'linewidth', 1.5 );
p(iSub+2) = bar(this_time, median(ko_mean)'*100, 'facecolor', 'k', 'facealpha',0.7);
p(iSub+1) = bar(this_time, median(wt_mean)'*100, 'facecolor', c_ord(3,:),'facealpha', 1);
ylabel('% NREM');
set(gca, 'xtick', 0:12:48, 'xticklabel', {'Z0', 'Z12', 'Z24', 'Z36', 'Z48'})
ylim([0 100]); xlim([0 48]);
legend([p(iSub+1) p(iSub+2)],{'wt', 'ko'}, 'box', 'off', 'orientation', 'horizontal');


subplot(3,4,9:11)
hold on
rectangle('position', [0 0.1 12 100], 'facecolor',[c_ord(5,:), 0.3], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [12 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.1], 'edgecolor', [.02 .02 .02 0]);
rectangle('position', [24 0.1 12 100], 'facecolor',[c_ord(5,:), 0.3], 'edgecolor', [c_ord(5,:), 0]);
rectangle('position', [36 0.1 12 100], 'facecolor',[0.2,0.2,0.2, 0.1], 'edgecolor', [.02 .02 .02 0]);
ko_mean = []; wt_mean = [];
for iSub = 1:length(Subjects)
    interval = mode(diff(Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.mins/60));  % get the minute interval from the time vector; 
    shift_val = ceil((Sub_scores.(PARAMS.Subjects.(Subjects{iSub}).genotype).(Subjects{iSub}).min_scores.clock_start/3600))*(60/interval); 

    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'ko')
        this_val =circshift(Sub_scores.ko.(Subjects{iSub}).min_scores.REM, shift_val);
        this_time = circshift(Sub_scores.ko.(Subjects{iSub}).min_scores.mins/3600, shift_val);
        ko_mean = [ko_mean; this_val];
    end
    
    if contains(PARAMS.Subjects.(Subjects{iSub}).genotype, 'wt')
        this_val =circshift(Sub_scores.wt.(Subjects{iSub}).min_scores.REM, shift_val);
        this_time = circshift(Sub_scores.wt.(Subjects{iSub}).min_scores.mins/3600, shift_val);
        wt_mean = [wt_mean; this_val];
    end
end
% p(iSub+1) = plot(this_time, median(wt_mean)'*100,'-*', 'color','k', 'MarkerEdgeColor','k', 'linewidth',1.5 );
% p(iSub+2) = plot(this_time, median(ko_mean)'*100,'-x', 'color',c_ord(2,:), 'MarkerEdgeColor',c_ord(2,:), 'linewidth', 1.5 );
p(iSub+2) = bar(this_time, median(ko_mean)'*100, 'facecolor', 'k', 'facealpha',0.7);
p(iSub+1) = bar(this_time, median(wt_mean)'*100, 'facecolor', c_ord(2,:),'facealpha', 1);
ylabel('% REM');
set(gca, 'xtick', 0:12:48, 'xticklabel', {'Z0', 'Z12', 'Z24', 'Z36', 'Z48'})
ylim([0 40]); xlim([0 48]);
legend([p(iSub+1) p(iSub+2)],{'wt', 'ko'}, 'box', 'off', 'orientation', 'horizontal');




% subplot(5,4,12)
% hold on
% rectangle('position', [0, 0, 3, 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
% rectangle('position', [3, 0, 3, 100], 'facecolor',[0.2,0.2,0.2, 0.2], 'edgecolor', [.02 .02 .02 0])
% 
% % data = [NHE6_tbl.REM(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Light')*100; NHE6_tbl.REM(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Light')*100;...
% % %     NaN(size(NHE6_tbl.Wake(NHE6_tbl.Genotype == 'ko'& NHE6_tbl.Light_Dark == 'Light'))); ...
% %     NHE6_tbl.REM(NHE6_tbl.Genotype == 'wt' & NHE6_tbl.Light_Dark == 'Dark')*100; NHE6_tbl.REM(NHE6_tbl.Genotype == 'ko' & NHE6_tbl.Light_Dark == 'Dark')*100]; 
% % % groups = [ones(1,length(data)), 2*ones(1,length(data)),  4*ones(1,length(data)), 5*ones(1,length(data))]; 
% % boxplot(data', 'Notch','on', 'Positions', [1,2,4,5])
% data = [reshape(circshift(wt_mean, ceil(this_time(1)),2), 1, numel(wt_mean)),reshape(circshift(ko_mean, ceil(this_time(1)),2), 1, numel(ko_mean))] ; 
% groups = [repmat([ones(1,12) 3*ones(1,12)],1,2),repmat([ones(1,12) 3*ones(1,12)],1,2), repmat([2*ones(1,12) 4*ones(1,12)],1,2),repmat([2*ones(1,12) 4*ones(1,12)],1,2)]; % groups ([wtL1-12 wtL25-36] [koL1-12 koL25-36] [wtD1-12 wtD25-36] [koD1-12 koD25-36]
% boxplot(data'*100, groups, 'Notch','on', 'Positions', [1,2,4,5])
% % groups = [ones(1,12),2*ones(1,12)];
% % boxplot(WAKE_out*100, groups);
% set(gca, 'xticklabel', {'w/t', 'k/o','w/t', 'k/o'})
% h = findobj(gca,'Tag','Box');
% patch(get(h(1),'XData'),get(h(1),'YData'),c_ord(2,:),'FaceAlpha',.5);
% patch(get(h(3),'XData'),get(h(3),'YData'),c_ord(2,:),'FaceAlpha',.5);
% patch(get(h(2),'XData'),get(h(2),'YData'),'k','FaceAlpha',.5);
% patch(get(h(4),'XData'),get(h(4),'YData'),'k','FaceAlpha',.5);
% ylim([0 20])
% 


% day_idx = nearest_idx(day_secs,all_tvec);
%
% day1_tvec = all_tvec(day_idx(1):day_idx(2));
% day2_tvec = all_tvec(day_idx(2)+1:day_idx(3));

% hour_idx = nearest_idx((5:24+5)*3600, day1_tvec)

% floor(max(all_tvec/3600))/24


