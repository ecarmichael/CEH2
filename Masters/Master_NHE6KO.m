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
PARAMS.Subjects.M05.LFP_Chan = 'CSCXX.ncs';  % best LFP channel.
PARAMS.Subjects.M05.EMG_Chan = 'CSCXX.ncs';  % best EMG channel.
PARAMS.Subjects.M05.genotype = 'wt';

% K/O
PARAMS.Subjects.M11.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M11.LFP_Chan = 'CSC39.ncs';  % best LFP channel.
PARAMS.Subjects.M11.EMG_Chan = 'CSC47.ncs';  % best EMG channel.
PARAMS.Subjects.M11.genotype = 'ko';

PARAMS.Subjects.M12.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M12.LFP_Chan = 'CSC1.ncs';  % best LFP channel.
PARAMS.Subjects.M12.EMG_Chan = 'CSC16.ncs';  % best EMG channel.
PARAMS.Subjects.M12.genotype = 'ko';

PARAMS.Subjects.M13.dir_name = 'NHE6KO_12_13_11_02'; % which foler type.  Can be either NHE6KO_12_13_11_02 or NHE6KO_05 based on which NLX system was used.
PARAMS.Subjects.M13.LFP_Chan = 'CSC17.ncs';  % best LFP channel.
PARAMS.Subjects.M13.EMG_Chan = 'CSC32.ncs';  % best EMG channel.
PARAMS.Subjects.M13.genotype = 'ko';
%% load and append data across multiple recording blocks.  To prevent buffer errors and corrupted data the 72hr recording period was broken down into several blocks of 8-11 hours (30hrs for mouse #05).


NH_preprocess(PARAMS.raw_data_dir, PARAMS.inter_dir); % run this script once to load, preprocess, and save all the 'good' channels for each subject to the PARAMS.inter_dir.


%% Classify the sleep states by visual inspection.
Subjects = fieldnames(PARAMS.Subjects);

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
        
        % catch cases where the sleep score will add in too many points to
        % match the range.
        if length(score{end}) ~= length(this_csc{iRec}.tvec(blocks(iB):end))
            %             score{end} =
        end
        
        
        % save the scores
        sleep_score{iRec}.tvec = this_csc{iRec}.tvec;
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
%% split the data into 1 hour blocks

% normalize the time across recording blocks
t_start = this_csc{1}.tvec(1);
Z0 = datetime(datestr('04-Sep-2020 08:00:00', 'dd-mmm-yyyy HH:MM:SS'));
clock_start = seconds(this_csc{1}.tstart - Z0);

all_tvec = [];
all_score = [];
all_data = []; 
Fs = mode(diff(sleep_score{iRec}.tvec));
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
    all_data = [all_data, [this_csc{iRec}.data(1,:); N_pad']]; 
    
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


%% break into 24 periods
% day_secs = (ceil(clock_start/3600))*3600:(24*3600):all_tvec(end);
hrs = (ceil(clock_start/3600))*3600:3600:52*3600;

for iH = 1:length(hrs)
    this_h_idx = nearest_idx([hrs(iH),hrs(iH)+3600],all_tvec); 
    disp(hrs(iH))
    WAKE(iH) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 1)/length(all_score(this_h_idx(1):this_h_idx(2)));
    NREM(iH) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 2)/length(all_score(this_h_idx(1):this_h_idx(2)));
    REM(iH) = sum(all_score(this_h_idx(1):this_h_idx(2)) == 3)/length(all_score(this_h_idx(1):this_h_idx(2)));
end

hour_scores.hrs = hrs; 
hour_scores.WAKE = WAKE;
hour_scores.NREM = NREM;
hour_scores.REM = REM; 

   save([PARAMS.inter_dir filesep Subjects{iSub} '_hour_data.mat'], 'hour_scores', '-v7.3')


% split into n x 24 array
hrs_out = circshift(hrs/3600,ceil(clock_start/3600));
WAKE_out = mean(reshape(circshift(WAKE,ceil(clock_start/3600)),24,2),2)';
NREM_out = mean(reshape(circshift(NREM,ceil(clock_start/3600)),24,2),2)';
REM_out = mean(reshape(circshift(REM,ceil(clock_start/3600)),24,2),2)';

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
set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
% set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
ylim([0 100]); xlim([1 24]);  

% bar plot for light vs dark
subplot(4,4,4)
groups = [ones(1,12),2*ones(1,12)];
h{1} = boxplot(WAKE_out, groups);
set(gca, 'xticklabel', {'light', 'dark'})

figure(1001)
subplot(4,4,5:7)
hold on
rectangle('position', [1 0.1 12 100], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(1:24, NREM_out*100,'-o', 'color', c_ord(3,:))
xlabel('time from light onset')
ylabel('% NREM');
% set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex');
ylim([0 100]); xlim([1 24]); 

% bar plot for light vs dark
subplot(4,4,8)
groups = [ones(1,12),2*ones(1,12)];
h{1} = boxplot(NREM_out, groups);
set(gca, 'xticklabel', {'light', 'dark'})


figure(1001)
subplot(4,4,9:11)
hold on
rectangle('position', [1 0.1 12 15], 'facecolor',[c_ord(5,:), 0.5], 'edgecolor', [c_ord(5,:), 0])
plot(1:24, REM_out*100,'-o', 'color', c_ord(2,:))
xlabel('time from light onset')
ylabel('% REM');
set(gca, 'xtick', 1:6:24, 'XTickLabel',XTickString,'TickLabelInterpreter','latex')
% set(gca, 'xtick', 1:6:24, 'xticklabel', {'Z0', 'Z6', 'Z12', 'Z18', 'Z24'})
% set(gca, 'xtick', 1:2:24, 'xticklabel', {'Z0','Z2', 'Z4', 'Z6','Z8', 'Z10', 'Z12','Z14', 'Z16','Z18','Z20', 'Z22', 'Z24'})
ylim([0 15]); xlim([1 24]); 

subplot(4,4,12)
groups = [ones(1,12),2*ones(1,12)];
h{1} = boxplot(REM_out, groups);
set(gca, 'xticklabel', {'light', 'dark'})


%% break out the sleep states and get the normalized PSD

%%% NEEDS all_data!! 
all_rem = all_data;
all_rem(all_score ~= 3) = []; 

win = 2048; 
[REM_ppx, REM_F] = pwelch(all_rem, hanning(win), win/2, 2*win, Fs);

subplot(4,4,13)
plot(REM_F, 10*log10(REM_ppx)); 
xlim([0 20])



% day_idx = nearest_idx(day_secs,all_tvec);
% 
% day1_tvec = all_tvec(day_idx(1):day_idx(2));
% day2_tvec = all_tvec(day_idx(2)+1:day_idx(3));

% hour_idx = nearest_idx((5:24+5)*3600, day1_tvec)

% floor(max(all_tvec/3600))/24


