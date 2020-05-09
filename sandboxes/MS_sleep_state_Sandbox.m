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

cfg_LFP = [];
cfg_LFP.fc = { 'CSC6.ncs'}; 
cfg_LFP.desired_sampling_frequency = 2000;

lfp = MS_LoadCSC(cfg_LFP);

cfg_EMG = [];
cfg_EMG.fc = { 'CSC1.ncs'}; 
cfg_EMG.desired_sampling_frequency = 2000;

emg = MS_LoadCSC(cfg_EMG);
%% filter
% delta filter.
cfg_filt_d = [];
cfg_filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_d.f  = [1 5];
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
% filter into the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [1 50];
cfg_filt_t.order = 12; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
wide_csc = FilterLFP(cfg_filt_t, csc);

%% collect the channels and save them for faster debugging. 
csc = emg;

csc.data(5,:) = wide_csc.data; 
csc.data(4,:) = theta_csc.data; 
csc.data(3,:) = delta_csc.data; 
csc.data(2,:) = lfp.data;

csc.label = {'emg', 'lfp', 'delta', 'theta'}; 
 clear theta_csc delta_csc emg lfp

 save('D:\Dropbox (Williams Lab)\Sleep_state\sleep_lfp_540day5.mat', 'csc'); 

%% load some premade lfp data to save time. 
 load('D:\Dropbox (Williams Lab)\Sleep_state\sleep_lfp_540day5.mat', 'csc'); 

csc_detect = csc;
csc_detect.data = [];
csc_detect.label = {};
 % convert emg to root mean square
csc_detect.data(1,:) = rms(csc.data(1,:)); 
csc_detect.label{1} = 'emg_rms'; 
 
 
 % normalize delta to 1-50Hz
 csc_detect.data(2,:) = csc.data(3,:) ./ bandpower(csc.data(2,:), csc.cfg.hdr{1}.SamplingFrequency, [1 50]); 
 
 csc_detect.data(3,:) = csc.data(3,:) ./ csc_
 
 %% bandpower method
 
 wide = bandpower(csc.data(3,:)); 
 
 %% short fourier transform
 f_d = [1 4];
 f_t = [6 9];
 f_w = [1 50]; 
 
%  [S, F, T] = stft(csc.data(2,:), csc.cfg.hdr{1}.SamplingFrequency, 'Window', hamming(csc.cfg.hdr{1}.SamplingFrequency*40));
 
 [~, F, T,P] =  spectrogram(csc.data(2,:),rectwin(csc.cfg.hdr{1}.SamplingFrequency*8),csc.cfg.hdr{1}.SamplingFrequency*4,0.5:.5:80,csc.cfg.hdr{1}.SamplingFrequency);
 figure
 subplot(4,1,1)
 imagesc(T,F,10*log10(P))
 axis xy
 
 subplot(4,1,2)
 plot(csc.tvec - csc.tvec(1), csc.data(1,:))
 
 
  subplot(4,1,3)
 plot(csc.tvec - csc.tvec(1), rms(csc.data(3,:)))
 
 
  subplot(4,1,4)
 plot(csc.tvec - csc.tvec(1), csc.data(4,:))
 
%  smag = mag2db(abs(S));
%         pcolor(seconds(T),F(find(F==0):find(F == 120)),smag(find(F==0):find(F == 120), :))
%         xlabel('Time (s)')
%         ylabel('Frequency (Hz)')
%         shading flat
%         colorbar
%         ylim([0 120])
  
        
        for iS = size(smag,2):-1:1
           dw_r(iS) = mean(smag(find(F == f_d(1)):find(F == f_d(2)),iS))/ mean(smag(find(F == f_w(1)):find(F == f_w(2)),iS));
           td_r(iS) = mean(smag(find(F == f_t(1)):find(F == f_t(2)),iS))/ mean(smag(find(F == f_d(1)):find(F == f_d(2)),iS));
           td_r(iS) = mean(smag(find(F == f_t(1)):find(F == f_t(2)),iS))/ mean(smag(find(F == f_d(1)):find(F == f_d(2)),iS));

            
            
        end
 
 