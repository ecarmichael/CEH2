%% MS_Sleep_detect_sandbox



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
    PARAMS.code_seqnmf_dir = 'C:\Users\ecarm\Documents\GitHub\seqNMF'; % where the multisite repo can be found
    
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


%% load some LFP data
cd('D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LFP\day5\2019-12-15_12-29-53_540day5base2');

cfg_LFP = [];
cfg_LFP.fc = { 'CSC6.ncs'};
cfg_LFP.desired_sampling_frequency = 2000;

lfp = MS_LoadCSC(cfg_LFP);


% cfg_LFP = [];
% cfg_LFP.fc = { 'CSC7.ncs'};
% cfg_LFP.desired_sampling_frequency = 2000;
% 
% lfp_7 = MS_LoadCSC(cfg_LFP);

cfg_EMG = [];
cfg_EMG.fc = { 'CSC1.ncs'};
cfg_EMG.desired_sampling_frequency = 2000;

emg = MS_LoadCSC(cfg_EMG);
% fix issue where first and last points are 0. makes smoothing an issue.


%% manually score the sleep states

cfg_sleep = [];
cfg_sleep.tvec_range = [0 10];  % number of seconds per window.
cfg_sleep.emg_range = [-0.001 0.001]; % default, should be based on real data.
cfg_sleep.emg_chan = 1; % emg channel.  Can be empty.
cfg_sleep.lfp_chans = 2; % lfp channels to be plotted can be empty. Can be 1 or more, but best to keep it less than 3. should be rows in csc.data.
cfg_sleep.state_name = {'Wake', 'SWS', 'REM', 'Quiescence', 'Transition', 'Redo'}; %

score = MS_Sleep_score_UI(cfg_sleep, csc.tvec, csc.data, abs(hilbert(emg.data(2,:))));



%% filter
% delta filter.
cfg_filt_d = [];
cfg_filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_d.f  = [1 4];
cfg_filt_d.order = 8; %type filter order
% cfg_filt_d.display_filter = 1; % use this to see the fvtool and wait for a keypress before continuing.
delta_csc = FilterLFP(cfg_filt_d,lfp);
close all


% filter into the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [6 9];
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
theta_csc = FilterLFP(cfg_filt_t, lfp);

% 'wide' for normalization
cfg_filt_t = [];
cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [1 50];
cfg_filt_t.order = 12; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
wide_csc = FilterLFP(cfg_filt_t, lfp);

% 'SWR' for SWS detection
cfg_filt_swr = [];
cfg_filt_swr.type = 'butter'; %Cheby1 is sharper than butter
cfg_filt_swr.f  = [140 250]; % broad, could use 150-200?
cfg_filt_swr.order = 4; %type filter order (fine for this f range)
% reason.  .
swr_csc = FilterLFP(cfg_filt_swr, lfp);

% bandpass the emg
cfg_filt_t = [];
cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [15 300];
cfg_filt_t.order = 16; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
emg = FilterLFP(cfg_filt_t, emg);


%% emg correction as per Buz lab code
% axx(1) =subplot(4,1,1);
% plot(emg.tvec, abs(hilbert(emg.data(1,:))))
% axx(2) =subplot(4,1,2);
%
% plot(emg.tvec, zscore(abs(hilbert(emg.data(1,:))).^2))

emg.data(2,:) = abs(hilbert(emg.data(1,:)));
emg.label{2} = 'Filt';
% axx(3) =subplot(4,1,3);
% plot(emg.tvec, emg.data(2,:));
% emg.data(3,:) = detrend(hilbert(emg.data(1,:)));
% axx(4) =subplot(4,1,4);
% plot(emg.tvec, emg.data(3,:));
% linkaxes(axx, 'x')
%
% %Min/Max Normalize
% EMG = bz_NormToRange(emg.data(2,:),[0 1]);
%
%
% dataspan = diff([min(emg.data(2,:)) max(emg.data(2,:))]);
% rangespan = diff([0 1]);
%
% normdata = (emg.data(2,:)-min(emg.data(2,:)))./dataspan;
% normdata = normdata.*rangespan+0;
%% collect the channels and save them for faster debugging.
csc = emg;
csc.data(7,:) = swr_csc.data;
csc.data(6,:) = wide_csc.data;
csc.data(5,:) = theta_csc.data;
csc.data(4,:) = delta_csc.data;
csc.data(3,:) = lfp.data;

emg_chan = 2;
raw_chan = 3; 
d_chan = 4;
t_chan = 5;
w_chan = 6;
r_chan = 7; 

csc.label = {'emg', 'emg_filt', 'lfp', 'delta', 'theta', 'wide', 'ripple'};
%  clear theta_csc delta_csc emg lfp

%  save('D:\Dropbox (Williams Lab)\Sleep_state\sleep_lfp_540day5.mat', 'csc');

%% load some premade lfp data to save time.
load('D:\Dropbox (Williams Lab)\Sleep_state\sleep_lfp_540day5.mat', 'csc');

% csc_detect = csc;
% csc_detect.data = [];
% csc_detect.label = {};
% convert emg to root mean square
% csc_detect.data(1,:) = rms(csc.data(1,:));
% csc_detect.label{1} = 'emg_rms';


% normalize delta to 1-50Hz
%  csc_detect.data(2,:) = csc.data(3,:) ./ bandpower(csc.data(2,:), csc.cfg.hdr{1}.SamplingFrequency, [1 50]);

%  csc_detect.data(3,:) = csc.data(3,:) ./ csc_

%  %% bandpower method
%
%  wide = bandpower(csc.data(3,:));


%% convert EMG into bins and get the rms
bins = 1:csc.cfg.hdr{1}.SamplingFrequency*2:length(csc.data(1,:));
emg_binned = NaN(size(bins));
%  delta_binned = NaN(size(bins));
%  theta_binned = NaN(size(bins));
%  wide_binned = NaN(size(bins));
%
%  delta_w_binned = delta_binned;
%  theta_w_binned = theta_binned;
%  wide_w_binned = wide_binned;
theta_binned = NaN(size(bins));
delta_binned =  NaN(size(bins));
emg_binned =  NaN(size(bins));
wide_binned = NaN(size(bins));
swr_binned =  NaN(size(bins));

t_d_binned = NaN(size(bins));
d_emg_binned =  NaN(size(bins));
t_emg_binned =  NaN(size(bins));
swr_t_binned = NaN(size(bins));
emg_d_binned = NaN(size(bins));
emg_t_binned = NaN(size(bins));
emg_swr_binned = NaN(size(bins));

tvec_binned = csc.tvec(bins) - csc.tvec(1); 

% if exist('score')
%     Sleepscore_binned = NaN(size(bins));
% end
tic;
for iB = length(bins):-1:1
    
    if iB == length(bins)
        
%         emg_binned(iB) = rms(emg.data(2,bins(1):end));
%         delta_binned(iB) = rms(csc.data(3,bins(1): end));
%         theta_binned(iB) = rms(csc.data(4,bins(1): end));
%         wide_binned(iB) = rms(csc.data(5,bins(1): end));
%         swr_binned(iB) = rms(csc.data(5,bins(1): end));
%         delta_w_binned(iB) = bandpower(csc.data(2,bins(iB): end), csc.cfg.hdr{1}.SamplingFrequency,[ 1,4]);
%         theta_w_binned(iB) = bandpower(csc.data(2,bins(iB): end), csc.cfg.hdr{1}.SamplingFrequency,[ 6,9]);
%         wide_w_binned(iB) = bandpower(csc.data(2,bins(iB): end), csc.cfg.hdr{1}.SamplingFrequency,[ 1,50]);
%         swr_w_binned(iB) =  bandpower(csc.data(2,bins(iB): end), csc.cfg.hdr{1}.SamplingFrequency,[ 125,250]);
        %


        d_w_binned(iB) = mean(abs(hilbert(csc.data(d_chan,bins(iB):end)))./abs(hilbert(csc.data(w_chan,bins(iB):end))));
        emg_binned(iB) = mean(csc.data(emg_chan,bins(iB):end));
        t_d_binned(iB) = mean(abs(hilbert(csc.data(t_chan,bins(iB):end)))./abs(hilbert(csc.data(d_chan,bins(iB):end))));
        d_emg_binned(iB) =  mean(abs(hilbert(csc.data(d_chan,bins(iB):end)))./csc.data(emg_chan,bins(iB):end));
        t_emg_binned(iB) =  mean(abs(hilbert(csc.data(t_chan,bins(iB):end)))./csc.data(emg_chan,bins(iB):end));
        t_swr_binned(iB) =  mean(abs(hilbert(csc.data(t_chan,bins(iB):end)))./abs(hilbert(csc.data(r_chan,bins(iB):end))));
        swr_t_binned(iB) =  mean(abs(hilbert(csc.data(r_chan,bins(iB):end)))./abs(hilbert(csc.data(t_chan,bins(iB):end))));
        swr_emg_binned(iB) =  mean(abs(hilbert(csc.data(r_chan,bins(iB):end)))./abs(hilbert(csc.data(emg_chan,bins(iB):end))));

        emg_d_binned(iB) = mean(csc.data(emg_chan,bins(iB):end)./abs(hilbert(csc.data(d_chan,bins(iB):end))));
        emg_t_binned(iB) = mean(csc.data(emg_chan,bins(iB):end)./abs(hilbert(csc.data(t_chan,bins(iB):end))));
        emg_swr_binned(iB) = mean(csc.data(emg_chan,bins(iB):end)./abs(hilbert(csc.data(r_chan,bins(iB):end))));
        
        t_d_emg_binned(iB) = mean((abs(hilbert(csc.data(t_chan,bins(iB):end)))./abs(hilbert(csc.data(emg_chan,bins(iB):end))))./...
            (abs(hilbert(csc.data(d_chan,bins(iB):end)))./abs(hilbert(csc.data(emg_chan,bins(iB):end)))));
                
%         t_d_emg_binned(iB) = t_d_binned(iB)./mean(csc.data(emg_chan,bins(iB):end));
%         
        
    else
        
%          emg_binned(iB) = rms(emg.data(2,bins(1): bins(iB+1)));
%         delta_binned(iB) = rms(csc.data(3,bins(1): bins(iB+1)));
%         theta_binned(iB) = rms(csc.data(4,bins(1): bins(iB+1)));
%         wide_binned(iB) = rms(csc.data(5,bins(1): bins(iB+1)));
%         swr_binned(iB) = rms(csc.data(5,bins(1): bins(iB+1)));
                
        
%         delta_w_binned(iB) = bandpower(csc.data(raw_chan,bins(iB): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 1,4]);
%         theta_w_binned(iB) = bandpower(csc.data(raw_chan,bins(iB): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 6,9]);
%         wide_w_binned(iB) = bandpower(csc.data(raw_chan,bins(iB): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 1,50]);
%         swr_w_binned(iB) =  bandpower(csc.data(raw_chan,bins(iB): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 125,250]);
        %
        d_w_binned(iB) = mean(abs(hilbert(csc.data(d_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(w_chan,bins(iB):bins(iB+1)))));
        emg_binned(iB) = mean(csc.data(emg_chan,bins(iB):bins(iB+1)));
        t_d_binned(iB) = mean(abs(hilbert(csc.data(t_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(d_chan,bins(iB):bins(iB+1)))));
        d_emg_binned(iB) =  mean(abs(hilbert(csc.data(d_chan,bins(iB):bins(iB+1))))./csc.data(emg_chan,bins(iB):bins(iB+1)));
        t_emg_binned(iB) =  mean(abs(hilbert(csc.data(t_chan,bins(iB):bins(iB+1))))./csc.data(emg_chan,bins(iB):bins(iB+1)));
        t_swr_binned(iB) =  mean(abs(hilbert(csc.data(t_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(r_chan,bins(iB):bins(iB+1)))));
        swr_t_binned(iB) =  mean(abs(hilbert(csc.data(r_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(t_chan,bins(iB):bins(iB+1)))));
        swr_emg_binned(iB) =  mean(abs(hilbert(csc.data(r_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(emg_chan,bins(iB):bins(iB+1)))));

        emg_d_binned(iB) = mean(csc.data(emg_chan,bins(iB):bins(iB+1))./abs(hilbert(csc.data(d_chan,bins(iB):bins(iB+1)))));
        emg_t_binned(iB) = mean(csc.data(emg_chan,bins(iB):bins(iB+1))./abs(hilbert(csc.data(t_chan,bins(iB):bins(iB+1)))));
        emg_swr_binned(iB) = mean(csc.data(emg_chan,bins(iB):bins(iB+1))./abs(hilbert(csc.data(r_chan,bins(iB):bins(iB+1)))));
        
        t_d_emg_binned(iB) = mean((abs(hilbert(csc.data(t_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(emg_chan,bins(iB):bins(iB+1)))))./...
            (abs(hilbert(csc.data(d_chan,bins(iB):bins(iB+1))))./abs(hilbert(csc.data(emg_chan,bins(iB):bins(iB+1))))));

    end
    %
    %
end
toc


% td_rms = theta_binned./delta_binned;
% d_emg_rms = delta_binned./emg_binned;
% t_emg_rms = theta_binned./emg_binned;
% swr_emg_rms = swr_binned./emg_binned;
% emg_d_rms = emg_binned./ delta_binned; 
% emg_t_rms = emg_binned./theta_binned;
% emg_swr_rms = emg_binned./theta_binned;



%  emg_d_binned = emg_d_binned(1:end-1);
%  emg_t_binned = emg_t_binned(1:end-1);
%  d_emg_binned = d_emg_binned(1:end-1);
%  t_emg_binned = t_emg_binned(1:end-1);
%  t_d_binned = t_d_binned(1:end-1);

%         emg_binned = emg_binned(1:end-2);
%           delta_binned = delta_binned(1:end-1);
%          theta_binned = theta_binned(1:end-1);
%          wide_binned = wide_binned(1:end-1);
%  delta_binned = 10*log10(delta_binned(2:end-1));
%  theta_binned = 10*log10(theta_binned(2:end-1));
%  wide_binned = 10*log10(wide_binned(2:end-1));
%
%  delta_w_binned = 10*log10(delta_w_binned(2:end-1));
%  theta_w_binned = 10*log10(theta_w_binned(2:end-1));
%  wide_w_binned = 10*log10(wide_w_binned(2:end-1));


% %% get the sleep state for the blocks?
%  B_SW_State = NaN(size(t_d_binned));
%  B_REM_State = NaN(size(t_d_binned));
%  B_W_State = NaN(size(t_d_binned));
%  B_S_State = NaN(size(t_d_binned));
%
%   B_W_State(emg_d_binned > .4 & t_emg_binned > 0.5) = 1; % get putative Wake
%  B_REM_State(t_emg > 7) = 1; % get putative REM
%  B_SW_State(d_emg > 6 & t_emg < 3) = 1;  % get putative SWS
%  S_State(isnan(S_State)) = 3; % get unclassified.

%% plot the bandpower versons
sample_multi = 2;

[~, F, T,P] =  spectrogram(csc.data(raw_chan,:),rectwin(csc.cfg.hdr{1}.SamplingFrequency*4),csc.cfg.hdr{1}.SamplingFrequency*sample_multi,0.5:.1:14,csc.cfg.hdr{1}.SamplingFrequency);

% to do: try normalizing the data within the band to make transitions more
% clear than raw amplitude.
d_amp = abs(hilbert(csc.data(d_chan,:)));
t_amp = abs(hilbert(csc.data(t_chan,:)));
r_amp = abs(hilbert(csc.data(r_chan,:)));
emg_amp = emg.data(2,:);


%  d_amp = MS_norm_range(d_amp, 0,1);
%  t_amp = MS_norm_range(t_amp, 0,1);
%  r_amp = MS_norm_range(r_amp, 0,1);
%  emg_amp = MS_norm_range(emg_amp, 0,1);

tvec_down = csc.tvec(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end)-csc.tvec(1);
score_down = score(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);


t_d = smooth(t_amp./d_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
t_r = smooth(t_amp./r_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
r_t = smooth(r_amp./t_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';

d_amp = smooth(d_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
t_amp = smooth(t_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
r_amp = smooth(r_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
emg_amp = smooth(emg_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';

d_emg =  (d_amp./emg_amp);
t_emg =  (t_amp./emg_amp);
t_d_emg = ((t_amp./emg_amp)./(d_amp./emg_amp));
r_emg = (r_amp./emg_amp);

%  d_emg =  smooth(d_amp./emg_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
%  t_emg =  smooth(t_amp./emg_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
%  t_d_emg = smooth(t_amp./emg,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';

emg_d = (emg_amp./d_amp);
emg_t = (emg_amp./t_amp);
emg_swr = (emg_amp./r_amp);

%  emg_d = smooth(emg_amp./d_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
%  emg_t = smooth(emg_amp./t_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';
%  emg_swr = smooth(emg_amp./r_amp,csc.cfg.hdr{1}.SamplingFrequency*sample_multi,'moving')';


% downsample
t_d = t_d(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
d_emg = d_emg(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
t_emg = t_emg(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
r_emg = r_emg(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
d_amp = d_amp(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
t_amp = t_amp(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
r_amp = r_amp(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
emg_amp = emg_amp(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);
t_d_emg = t_d_emg(1:csc.cfg.hdr{1}.SamplingFrequency*sample_multi:end);

SW_State = NaN(size(score));
REM_State = NaN(size(score));
W_State = NaN(size(score));
Q_State = NaN(size(score));
T_State = NaN(size(score));

W_State(score ==1) = 1;
SW_State(score ==2) = 1;
REM_State(score ==3) = 1;
Q_State(score==4) = 1;
T_State(score==5) = 1;

%% try PCA of d t emg????
% [p_coeff, p_score] = pca([zscore(d_amp); zscore(t_amp);zscore(r_amp); zscore(emg_amp)]);

[p_coeff, p_score] = pca([zscore(d_emg); zscore(t_emg);zscore(r_amp); zscore(emg_amp)]);
figure(200)
subplot(411)
plot(T, emg_binned(1:end-2))
subplot(412)
plot(T, MS_norm_range(p_coeff(1:end-2,1), 0, 1))
subplot(413)
plot(T, p_coeff(1:end-2,2))
subplot(414)
plot(T, p_coeff(1:end-2,3))

figure(201)
scatter3(p_coeff(:,3),p_coeff(:,2), p_coeff(:,1), '.')
%% try some K means clustering
data_in = [d_emg_binned', t_emg_binned', t_d_binned',emg_binned'];

MS_kmean_scatter(data_in, 3,[4,3,2,1], 20);


%% try auto-detect   scores {'Wake', 'SWS', 'REM', 'Quiescence', 'Transition'}

SW_State_auto = NaN(size(t_emg));
REM_State_auto = NaN(size(t_emg));
W_State_auto = NaN(size(t_emg));
Q_State_auto = NaN(size(t_emg));
T_State_auto = NaN(size(t_emg));
All_State_auto = NaN(size(t_emg));

d_emg_th = 3; % if greater should be awake
t_emg_th = 5; % if greater, then REM
emg_t_th = 0.3;
emg_amp_th = 0.00003; % if greater Awake
emg_amp_th_low = 0.00002; % lower end. if less than SWS or REM. 
t_d_th = 1;
r_t_th = 2;
t_r_th = 0.6;
t_d_emg = 30000;

%  for td = 0:.1:7
% Rem (rare
All_State_auto((t_emg >= t_emg_th)  & (t_d > t_d_th) ) = 3;
% SWS events
All_State_auto((d_emg >= d_emg_th)  & (t_d > t_d_th) ) = 2;


figure(1010)
subplot(2,1,1)
hold on
%  plot(tvec_down, (emg_amp < emg_amp_th)*1, '.r');
plot(tvec_down, (d_emg < d_emg_th)*1, '.r');text(-1000, 1,'d emg');
plot(tvec_down, ((All_State_auto==3)*3)+.5, '.m')
plot(tvec_down, (t_emg >= t_emg_th)*2, '.g');  text(-1000, 2,'t emg');
plot(tvec_down, (t_d > t_d_th)*3, '.b');text(-1000,3,'t d');
% plot(tvec_down, (t_emg >= t_emg_th  & t_d > t_d_th)*4, '.k');text(-500, 4,'rem?');
plot(tvec_down, ((score_down==3)*3)+.25, '.k');

ylim([.8 5])
xlim([min(tvec_down) max(tvec_down)])

subplot(2,1,2)
hold on
plot(tvec_down, ((All_State_auto==3)*3)+.2, '.r');
plot(tvec_down, score_down, '.k');
text(-1000, 1, 'Wake');
text(-1000, 2, 'SWS');
text(-1000, 3, 'REM');
text(-1000, 4, 'Q');
text(-1000, 5, 'T');

xlim([min(csc.tvec-csc.tvec(1)) max(csc.tvec-csc.tvec(1))])

% check the score.
score_rem = score_down==3;

rem_auto = t_emg >= t_emg_th  & emg_amp < emg_amp_th & t_d > t_d_th;
overlap = nansum(rem_auto == score_rem')/length(score_rem);

title(['Overlap: ' num2str(overlap*100,4) '%'])

%  end
% SWS
%  All_State_auto(d_emg >= 15 & emg_amp_n < emg_amp_th & t_d < t_d_th) = 2;
%   % Quiescence
%   All_State_auto(emg_amp < emg_amp_th  & t_emg < t_emg_th) = 4;
%  % Wake
%   All_State_auto(emg_amp >= emg_amp_th  & t_d >= t_d_th) = 1;
%
% figure(222)
% hold on
% plot(csc.tvec-csc.tvec(1), score, '.k');
% plot(csc.tvec-csc.tvec(1), All_State_auto+0.2, '.r');
% ylim([0 6]);
% set(gca, 'ytick', 1:5, 'yticklabel', {'Wake', 'SWS', 'REM', 'Quiet', 'Tran'});
% %  % to do loop this to optimize.
% % REM_State_auto(t_emg > 5 & emg_amp_n < 0.2 ) = 1; % get putative REM
% % All_State_auto(~isnan(REM_State_auto)) =3;
%
% % SW_State_auto(d_emg > 15 & emg_amp_n < 0.2 & t_d < 2) = 1;  % get putative SWS
% % W_State_auto(emg_amp_n >= 0.2 & REM_State_auto ~=1) = 1; % get putative Wake
% % % Q_State_auto(isnan(S_State)) = 3; % get unclassified.
% overlap = nansum(All_State_auto == score)/length(score);
%
% title(['Overlap: ' num2str(overlap*100,2) '%'])


%% plot the power of the filtered signals.
c_ord = linspecer(5);

Subs = 9;
figure(115)
ax2(1) =subplot(Subs,1,1:2);
P_p = 10*log10(P);
imagesc(T,F,P_p)
axis xy
hold on
plot(csc.tvec-csc.tvec(1), W_State+max(F)+1,  'color', c_ord(1,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), Q_State+max(F)+1,  'color', c_ord(5,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), SW_State+max(F)+1, 'color', c_ord(3,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), REM_State+max(F)+1, 'color', c_ord(2,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), T_State+max(F)+1, 'color', c_ord(4,:), 'linewidth', 4)
ylim([F(1) F(end)+5])
ax2(2) =subplot(Subs,1,3);
plot(tvec_down, emg_amp)
text(-500, nanmedian(emg_amp),'emg/amp');
hline(emg_amp_th)

ax2(3) =subplot(Subs,1,4);
if length(d_amp) ~= length(csc.tvec)
    plot(tvec_down, d_amp)
else
    plot(csc.tvec-csc.tvec(1), d_amp)
end
text(-500, nanmedian(d_amp),'d/amp')
hline(d_amp_th)


ax2(4) =subplot(Subs,1,5);
if length(t_amp) ~= length(csc.tvec)
    plot(tvec_down, t_amp)
else
    plot(csc.tvec-csc.tvec(1), t_amp)
end
text(-500, nanmedian(t_amp),'t/amp')

ax2(5) =subplot(Subs,1,6);
if length(t_emg) ~= length(csc.tvec)
    plot(tvec_down, t_emg)
else
    plot(csc.tvec-csc.tvec(1), t_emg)
end
text(-500, nanmedian(t_emg),'t/emg')
hline(t_emg_th)


ax2(6) =subplot(Subs,1,7);
if length(d_emg) ~= length(csc.tvec)
    plot(tvec_down, d_emg)
else
    plot(csc.tvec-csc.tvec(1), d_emg)
end
text(-500, nanmedian(d_emg),'d/emg')
hline(d_emg_th)

ax2(7) =subplot(Subs,1,8);
if length(t_d) ~= length(csc.tvec)
    plot(tvec_down, t_d)
else
    plot(csc.tvec-csc.tvec(1), t_d)
end
text(-500, nanmedian(t_d),'t/d')
hline(t_d_th)


ax2(8) =subplot(Subs,1,9);
%  swr_t([1 end]) = NaN;
%  t_swr([1 end]) = NaN;
%  emg_d(end-10:end) = NaN;
if length(r_t) ~= length(csc.tvec)
    plot(tvec_down, r_t)
else
    plot(csc.tvec-csc.tvec(1), r_t)
end
text(-500, nanmedian(r_t),'r t')
hline(t_emg_th)
linkaxes(ax2, 'x');
xlim([T(1) T(end)])

%% try to auto classify. 


SW_State_auto_B = NaN(size(t_d_binned));
REM_State_auto_B = NaN(size(t_d_binned));
W_State_auto_B = NaN(size(t_d_binned));
Q_State_auto_B = NaN(size(t_d_binned));
T_State_auto_B = NaN(size(t_d_binned));
All_State_auto_B = NaN(size(t_d_binned));

d_emg_th = 6; % if greater should be awake
t_emg_th = 6; % if greater, then REM
emg_t_th = 0.3;
emg_amp_th = 0.00003; % if greater Awake
emg_amp_th_low = 0.00002; % lower end. if less than SWS or REM. 

t_d_th = 1.25;
r_t_th = 2;
t_r_th = 0.6;
t_d_emg = 30000;

%  for td = 0:.1:7
% Rem (rare
All_State_auto_B((t_emg_binned >= t_emg_th)  & (t_d_binned > t_d_th) & (emg_binned < emg_amp_th) ) = 3;
% SWS events
All_State_auto_B((d_emg_binned >= d_emg_th)  & (t_d_binned < t_d_th) ) = 2;
% Quiet
All_State_auto_B((d_emg_binned >= d_emg_th)  & (t_d_binned < t_d_th) ) = 4;
% Awake
All_State_auto_B((d_emg_binned < d_emg_th) & (t_emg_binned < t_emg_th)  & (emg_binned > emg_amp_th) ) = 1;


% close(1010); 
figure(1010)
subplot(2,1,1)
hold on
% REM tests
% %  plot(tvec_down, (emg_amp < emg_amp_th)*1, '.r');
% plot(tvec_binned, (emg_binned < emg_amp_th)*1, '.r');text(-1000, 1,' emg amp');
% plot(tvec_binned, ((All_State_auto_B==3)*3)+.5, '.m')
% plot(tvec_binned, (t_emg_binned >= t_emg_th)*2, '.g');  text(-1000, 2,'t emg');
% plot(tvec_binned, (t_d_binned > t_d_th)*3, '.b');text(-1000,3,'t d');
% plot(tvec_binned, (d_emg_binned < d_emg_th)*4, '.b');text(-1000,4,'d emg');
% plot(tvec_binned, ((score_down==3)*3)+.25, '.k');

% SWS
plot(tvec_binned, (t_d_binned < t_d_th)*2, '.g');  text(-1000, 2,'t emg');
plot(tvec_binned, (d_emg_binned >= d_emg_th)*4, '.b');text(-1000,4,'d emg');

plot(tvec_binned, ((score_down==2)*3)+.25, '.k');
plot(tvec_binned, ((All_State_auto_B==2)*3)+.5, '.m')

% awake
% plot(tvec_binned, (emg_binned > emg_amp_th)*1, '.r');text(-1000, 1,' emg amp');
% plot(tvec_binned, (t_emg_binned < t_emg_th)*2, '.g');  text(-1000, 2,'t emg');
% plot(tvec_binned, (d_emg_binned < d_emg_th)*4, '.b');text(-1000,4,'d emg');
% 
% plot(tvec_binned, ((score_down==1)*3)+.25, '.k');
% plot(tvec_binned, ((All_State_auto_B==1)*3)+.5, '.m')

% plot(tvec_down, (t_emg >= t_emg_th  & t_d > t_d_th)*4, '.k');text(-500, 4,'rem?');
ylim([0.8 5])
xlim([min(tvec_binned) max(tvec_binned)])

subplot(2,1,2)
hold on
plot(tvec_binned, All_State_auto_B+.25, '.r');
plot(tvec_binned, score_down, '.k');
text(-1000, 1, 'Wake');
text(-1000, 2, 'SWS');
text(-1000, 3, 'REM');
text(-1000, 4, 'Q');
text(-1000, 5, 'T');
ylim([0.5 5.5])

xlim([min(csc.tvec-csc.tvec(1)) max(csc.tvec-csc.tvec(1))])

% check the score.
score_rem = score_down==3;

rem_auto = t_emg >= t_emg_th  & emg_amp < emg_amp_th & t_d > t_d_th;
overlap = nansum(rem_auto == score_rem')/length(score_rem);

title(['Overlap: ' num2str(overlap*100,4) '%'])


%% binned data blocks.
c_ord = linspecer(5);
figure(120)

ax(1) =subplot(9,1,1:2);
P_p = 10*log10(P);
imagesc(T,F,P_p)
axis xy
hold on
plot(csc.tvec-csc.tvec(1), W_State+max(F),  'color', c_ord(1,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), Q_State+max(F),  'color', c_ord(5,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), SW_State+max(F), 'color', c_ord(3,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), REM_State+max(F), 'color', c_ord(2,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), T_State+max(F), 'color', c_ord(4,:), 'linewidth', 4)
ylim([F(1) F(end)+1.5])

ax(9) =subplot(9,1,3);
plot(tvec_binned, emg_binned)
text(-400, mean(emg_binned),'emg amp')
hline(emg_amp_th)

% ax(2) =subplot(9,1,4);
% plot(T, p_coeff(1:end-2,:));
% text(-400, mean(p_coeff(1:end-2,:), 'all'),'pca')
% legend({'pc1', 'pc2', 'pc3'},'location', 'northwest','Orientation','horizontal')
% hline(.5)
ax(2) =subplot(9,1,4);
plot(tvec_binned, d_swr_emg_binned);
text(-400, mean(d_swr_emg_binned),'d swr emg')
% hline(t_d__th)


ax(4) =subplot(9,1,5);
plot(tvec_binned, d_emg_binned);
text(-400, mean(d_emg_binned),'d/emg')
hline(d_emg_th)

ax(5) =subplot(9,1,6);
plot(tvec_binned, t_d_binned)%, T, p_coeff(1:end-2,2)*100);
text(-400, mean(t_d_binned),'th /d')
% legend({'t/d','pc2'},'location', 'northwest','Orientation','horizontal')
hline(t_d_th)

ax(6) =subplot(9,1,7);
plot(tvec_binned, t_emg_binned)%, T, p_coeff(1:end-2,2)*200);
text(-400, mean(t_emg_binned),'th /emg')
% legend({'t/emg', 'pc2'},'location', 'northwest','Orientation','horizontal')
hline(t_emg_th)

ax(7) =subplot(9,1,8);
plot(tvec_binned, emg_t_binned);
text(-400, mean(emg_t_binned),'emg /t')
hline(emg_t_th)

ax(8) =subplot(9,1,9);
plot(tvec_binned, t_d_wide_binned);
text(-400, mean( t_d_wide_binned),'t/emg ./ d/emg')
hline(swr_t_th)

linkaxes(ax, 'x')
xlim([T(1) T(end)])


%% RMS binned data blocks.
c_ord = linspecer(5);
figure(122)

ax(1) =subplot(9,1,1:2);
P_p = 10*log10(P);
imagesc(T,F,P_p)
axis xy
hold on
plot(csc.tvec-csc.tvec(1), W_State+max(F),  'color', c_ord(1,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), Q_State+max(F),  'color', c_ord(5,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), SW_State+max(F), 'color', c_ord(3,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), REM_State+max(F), 'color', c_ord(2,:), 'linewidth', 4)
plot(csc.tvec-csc.tvec(1), T_State+max(F), 'color', c_ord(4,:), 'linewidth', 4)
ylim([F(1) F(end)+1.5])

ax(9) =subplot(9,1,3);
plot(tvec_binned, emg_binned)
text(-400, mean(emg_binned),'emg amp')
hline(emg_amp_th)

% ax(2) =subplot(9,1,4);
% plot(T, p_coeff(1:end-2,:));
% text(-400, mean(p_coeff(1:end-2,:), 'all'),'pca')
% legend({'pc1', 'pc2', 'pc3'},'location', 'northwest','Orientation','horizontal')
% hline(.5)
ax(2) =subplot(9,1,4);
plot(tvec_binned, d_emg_rms);
text(-400, mean(d_emg_rms),'d emg rms')
% hline(t_d__th)


ax(4) =subplot(9,1,5);
plot(tvec_binned, t_emg_rms);
text(-400, mean(t_emg_rms),'t_emg_rms')
hline(d_emg_th)

ax(5) =subplot(9,1,6);
plot(tvec_binned, t_d_)%, T, p_coeff(1:end-2,2)*100);
text(-400, mean(t_d_binned),'th /d')
% legend({'t/d','pc2'},'location', 'northwest','Orientation','horizontal')
hline(t_d_th)

ax(6) =subplot(9,1,7);
plot(tvec_binned, t_emg_binned)%, T, p_coeff(1:end-2,2)*200);
text(-400, mean(t_emg_binned),'th /emg')
% legend({'t/emg', 'pc2'},'location', 'northwest','Orientation','horizontal')
hline(t_emg_th)

ax(7) =subplot(9,1,8);
plot(tvec_binned, emg_t_binned);
text(-400, mean(emg_t_binned),'emg /t')
hline(emg_t_th)

ax(8) =subplot(9,1,9);
plot(tvec_binned, t_d_wide_binned);
text(-400, mean( t_d_wide_binned),'t/emg ./ d/emg')
hline(swr_t_th)

linkaxes(ax, 'x')
xlim([T(1) T(end)])


%% %%%%%%%%%%%%%%%%%%%%%%%
% short fourier transform method
f_d = [20 80];
f_t = [5 9];
f_w = [10 50];

%  [S, F, T] = stft(csc.data(2,:), csc.cfg.hdr{1}.SamplingFrequency, 'Window', hamming(csc.cfg.hdr{1}.SamplingFrequency*40));

[~, F, T,P] =  spectrogram(csc.data(raw_chan,:),hanning(csc.cfg.hdr{1}.SamplingFrequency*4),csc.cfg.hdr{1}.SamplingFrequency*2,1:1:150,csc.cfg.hdr{1}.SamplingFrequency);
%%
P_p = 10*log10(P);
% P_p = P_p.^2;
P_p = P_p./min(min(P_p));
figure(125)
ax(1) = subplot(9,1,1:4);
imagesc(T,F,P_p)
axis xy
d_p = []; w_p = []; t_p=[];

for iS = size(P_p,2):-1:1
    d_p(iS) = mean(P_p(find(F == f_d(1)):find(F == f_d(2)),iS));
    w_p(iS) = mean(P_p(find(F == f_w(1)):find(F == f_w(2)),iS));
    t_p(iS) = mean(P_p(find(F == f_t(1)):find(F == f_t(2)),iS));
end
d_norm = d_p./w_p;
t_norm = t_p./w_p;

td_norm = t_norm./d_norm;
td = t_p./d_p; 

%  subplot(5,1,2)
%  plot(csc.tvec - csc.tvec(1), csc.data(1,:))
%
ax(2) = subplot(9,1,5);
norm_t = csc.tvec - csc.tvec(1);
plot(norm_t(bins(2:end-1)),emg_binned(2:end-1))
xlabel('EMG rms')

ax(3) = subplot(9,1,6);
plot(T, d_p);
xlabel('1-4')

ax(4) =   subplot(9,1,7);
plot(T, d_norm);
xlabel('1-4/ 1-50')

ax(6) = subplot(9,1,8);
plot(T, td);
xlabel('6-9 / 1-4')

ax(7) = subplot(9,1,9);
plot(T, td_norm);
xlabel('6-9 / 1-4')
linkaxes(ax, 'x');

%   subplot(5,1,4)
%  plot(csc.tvec - csc.tvec(1), csc.data(3,:))


%   subplot(5,1,5)
%  plot(csc.tvec - csc.tvec(1), csc.data(4,:))
%
%  smag = mag2db(abs(S));
%         pcolor(seconds(T),F(find(F==0):find(F == 120)),smag(find(F==0):find(F == 120), :))
%         xlabel('Time (s)')
%         ylabel('Frequency (Hz)')
%         shading flat
%         colorbar
%         ylim([0 120])
%% check overlap. try model classification ??

% class_data(:,1) = score_down;
% class_data(:,2) = t_d';
% class_data(:,3) = emg_amp';
% class_data(:,4) = d_amp';
% class_data(:,5) = t_amp';
% class_data(:,6) = d_emg';
% class_data(:,7) = t_emg';
% class_data(:,8) = t_d_emg';

% best case was 77% with Medium Tree or linear discrim.  Most model types were around 75%.
%% try to auto class based on binned data
% spec_class= cat(2,score_down(2:end-1),emg_binned(2:end-1)', P_p');

% got to 75% with linear SVM  good at REM.  


%% try something that has already been scored.
cfg = [];
cfg.fc = {'J:\Williams_Lab\EV_testcalcium\CSC4.ncs'};
cfg.desired_sampling_frequency = 2000;
csc = MS_LoadCSC(cfg);

% load hypnogram
load('Hypnogram.mat');

%%
figure(11)
subplot(2,1,1)

plot(csc.tvec - csc.tvec(1), csc.data)
subplot(2,1,2)

plot(1:length(Hypnogram),Hypnogram)


%% morlet wavelet


% Morlet wavelet
Freq = 0.5:0.5:40;
coefsi = cwt(csc.data ,centfrq('cmor1 -1 ')*csc.cfg.hdr{1}.SamplingFrequency./Freq,'cmor1 -1 ');
figure;
subplot(2,1,1)
wm_image = imagesc(csc.tvec, Freq, abs(coefsi));

axis xy

caxis([ 0 0.0005])


[~, F, T,P] =  spectrogram(csc.data(1,:),rectwin(csc.cfg.hdr{1}.SamplingFrequency*8),csc.cfg.hdr{1}.SamplingFrequency*4,0.5:.5:40,csc.cfg.hdr{1}.SamplingFrequency);

subplot(2,1,2)
imagesc(T,F,10*log10(P))
axis xy


%% try SeqNMF [doesn't work in most basic form. only played with sampling (f and t)]
addpath(genpath(PARAMS.code_seqnmf_dir));
sample_multi = 4;
[~, F, T,P] =  spectrogram(csc.data(2,:),rectwin(csc.cfg.hdr{1}.SamplingFrequency*(2*sample_multi)),csc.cfg.hdr{1}.SamplingFrequency*sample_multi,0.5:.5:80,csc.cfg.hdr{1}.SamplingFrequency);


data_in = 10*log10(P)./max(10*log10(P));
Fs = csc.cfg.hdr{1}.SamplingFrequency/sample_multi;
splitN = floor(size(data_in,2)*.90);
trainNEURAL = data_in(:,1:splitN);
testNEURAL = data_in(:,(splitN+1):end);

rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;

K = 5;
L = 1; % units of seconds
Lneural = ceil(L*Fs);
% Lsong = ceil(L*SONGfs);
shg
display('Running seqNMF on real neural data (from songbird HVC, recorded by Emily Mackevicius, Fee Lab)')
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
    'lambdaL1W', .1, 'lambda', .005, 'maxiter', 100, 'showPlot', 1,...
    'lambdaOrthoW', 0);

p = .05; % desired p value for factors

display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(testNEURAL,W,p);

W = W(:,is_significant,:);
H = H(is_significant,:);
