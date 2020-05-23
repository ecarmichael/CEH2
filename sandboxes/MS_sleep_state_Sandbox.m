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


cfg_LFP = [];
cfg_LFP.fc = { 'CSC7.ncs'}; 
cfg_LFP.desired_sampling_frequency = 2000;

lfp_7 = MS_LoadCSC(cfg_LFP);

cfg_EMG = [];
cfg_EMG.fc = { 'CSC1.ncs'}; 
cfg_EMG.desired_sampling_frequency = 2000;

emg = MS_LoadCSC(cfg_EMG);
% fix issue where first and last points are 0. makes smoothing an issue. 
emg.data(1,1) = emg.data(1,2);
emg.data(1,end) = emg.data(1,end-1); 

%% emg correction as per Buz lab code
axx(1) =subplot(4,1,1);
plot(emg.tvec, abs(hilbert(emg.data(1,:))))
axx(2) =subplot(4,1,2);

plot(emg.tvec, zscore(abs(hilbert(emg.data(1,:))).^2))

emg.data(2,:) = smooth(abs(hilbert(emg.data(1,:))),emg.cfg.hdr{1}.SamplingFrequency*4,'moving');
axx(3) =subplot(4,1,3);
plot(emg.tvec, emg.data(2,:));
emg.data(3,:) = detrend(hilbert(emg.data(1,:)));
axx(4) =subplot(4,1,4);
plot(emg.tvec, emg.data(3,:));
linkaxes(axx, 'x')

%Min/Max Normalize
EMG = bz_NormToRange(emg.data(2,:),[0 1]);


dataspan = diff([min(emg.data(2,:)) max(emg.data(2,:))]);
rangespan = diff([0 1]);

normdata = (emg.data(2,:)-min(emg.data(2,:)))./dataspan;
normdata = normdata.*rangespan+0;
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
wide_csc = FilterLFP(cfg_filt_t, lfp);


% bandpass the emg
cfg_filt_t = [];
cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [15 300];
cfg_filt_t.order = 32; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .
emg = FilterLFP(cfg_filt_t, emg);
%% collect the channels and save them for faster debugging. 
csc = emg;

csc.data(5,:) = wide_csc.data; 
csc.data(4,:) = theta_csc.data; 
csc.data(3,:) = delta_csc.data; 
csc.data(2,:) = lfp.data;

csc.label = {'emg', 'lfp', 'delta', 'theta', 'wide'}; 
%  clear theta_csc delta_csc emg lfp

 save('D:\Dropbox (Williams Lab)\Sleep_state\sleep_lfp_540day5.mat', 'csc'); 

%% load some premade lfp data to save time. 
 load('D:\Dropbox (Williams Lab)\Sleep_state\sleep_lfp_540day5.mat', 'csc'); 

csc_detect = csc;
csc_detect.data = [];
csc_detect.label = {};
 % convert emg to root mean square
% csc_detect.data(1,:) = rms(csc.data(1,:)); 
% csc_detect.label{1} = 'emg_rms'; 
 
 
 % normalize delta to 1-50Hz
%  csc_detect.data(2,:) = csc.data(3,:) ./ bandpower(csc.data(2,:), csc.cfg.hdr{1}.SamplingFrequency, [1 50]); 
 
%  csc_detect.data(3,:) = csc.data(3,:) ./ csc_
 
 %% bandpower method
 
 wide = bandpower(csc.data(3,:)); 
 
 
 %% convert EMG into bins and get the rms
 bins = 1:csc.cfg.hdr{1}.SamplingFrequency*4:length(csc.data(1,:)); 
 emg_binned = NaN(size(bins)); 
%  delta_binned = NaN(size(bins));
%  theta_binned = NaN(size(bins));
%  wide_binned = NaN(size(bins));
%  
%  delta_w_binned = delta_binned;
%  theta_w_binned = theta_binned;
%  wide_w_binned = wide_binned;
 
 
 tic;
for iB = length(bins)-1:-1:1
        
%          emg_binned(iB) = rms(emg.data(2,bins(1): bins(iB+1)));
         delta_binned(iB) = rms(csc.data(3,bins(1): bins(iB+1))); 
         theta_binned(iB) = rms(csc.data(4,bins(1): bins(iB+1))); 
         wide_binned(iB) = rms(csc.data(5,bins(1): bins(iB+1))); 
         
%          delta_w_binned(iB) = bandpower(csc.data(2,bins(1): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 1,4]);
%          theta_w_binned(iB) = bandpower(csc.data(2,bins(1): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 6,9]);
%          wide_w_binned(iB) = bandpower(csc.data(2,bins(1): bins(iB+1)), csc.cfg.hdr{1}.SamplingFrequency,[ 1,50]);
% 
%      
%      
 end
 toc
 
        emg_binned = emg_binned(1:end-2);
          delta_binned = delta_binned(1:end-1); 
         theta_binned = theta_binned(1:end-1);
         wide_binned = wide_binned(1:end-1);
%  delta_binned = 10*log10(delta_binned(2:end-1));
%  theta_binned = 10*log10(theta_binned(2:end-1));
%  wide_binned = 10*log10(wide_binned(2:end-1));
%  
%  delta_w_binned = 10*log10(delta_w_binned(2:end-1));
%  theta_w_binned = 10*log10(theta_w_binned(2:end-1));
%  wide_w_binned = 10*log10(wide_w_binned(2:end-1));
 
 %% plot the bandpower versons
 

 [~, F, T,P] =  spectrogram(csc.data(2,:),rectwin(csc.cfg.hdr{1}.SamplingFrequency*8),csc.cfg.hdr{1}.SamplingFrequency*4,0.5:.5:50,csc.cfg.hdr{1}.SamplingFrequency);
 
 t_d = smooth(abs(hilbert(csc.data(4,:)))./abs(hilbert(csc.data(3,:))),csc.cfg.hdr{1}.SamplingFrequency*4,'moving');
 d_emg =  smooth(abs(hilbert(csc.data(3,:)))./emg.data(2,:),csc.cfg.hdr{1}.SamplingFrequency*4,'moving');
 t_emg =  smooth(abs(hilbert(csc.data(4,:)))./emg.data(2,:),csc.cfg.hdr{1}.SamplingFrequency*4,'moving');
 emg_d = smooth(emg.data(2,:)./abs(hilbert(csc.data(3,:))),csc.cfg.hdr{1}.SamplingFrequency*4,'moving');
 SW_State = NaN(size(t_d)); 
 REM_State = NaN(size(t_d)); 
 W_State = NaN(size(t_d)); 

  REM_State(t_emg > 5) = 1; % get putative REM 
 SW_State(d_emg > 6 & t_emg < 3) = 1;  % get putative SWS
 W_State(emg_d > .4 & REM_State ~=1) = 1; % get putative Wake 
 S_State(isnan(S_State)) = 3; % get unclassified. 
%  S_State(

 %% plot the power of the filtered signals. 
 c_ord = linspecer(4);
 
 Subs = 8;
 figure(115)
 ax2(1) =subplot(Subs,1,1:2);
 P_p = 10*log10(P);
 imagesc(T,F,P_p)
 axis xy
 hold on
 plot(csc.tvec-csc.tvec(1), SW_State+38,'*','markersize', 4, 'color', c_ord(4,:), 'linewidth', 4)
  plot(csc.tvec-csc.tvec(1), W_State+44,'*','markersize', 4, 'color', c_ord(1,:), 'linewidth', 4)
 plot(csc.tvec-csc.tvec(1), REM_State+41,'*','markersize', 4, 'color', c_ord(2,:), 'linewidth', 4)

 ax2(2) =subplot(Subs,1,3);
 plot(csc.tvec-csc.tvec(1), emg.data(2,:))
 xlabel('emg')
 ax2(3) =subplot(Subs,1,4);
 plot(csc.tvec-csc.tvec(1), smooth(abs(hilbert(csc.data(3,:))),csc.cfg.hdr{1}.SamplingFrequency*4,'moving'))
 xlabel('delta')
 
  ax2(4) =subplot(Subs,1,5);
 plot(csc.tvec-csc.tvec(1), emg_d)
 xlabel('delta/emg')
 
 ax2(6) =subplot(Subs,1,6);
 plot(csc.tvec-csc.tvec(1), smooth(abs(hilbert(csc.data(4,:))),csc.cfg.hdr{1}.SamplingFrequency*4,'moving'))
 xlabel('theta')
 
 ax2(7) =subplot(Subs,1,7);
 plot(csc.tvec-csc.tvec(1),t_d)
 xlabel('theta/delta')
 
 ax2(8) =subplot(Subs,1,8);
 plot(csc.tvec-csc.tvec(1),t_emg)
 xlabel('theta/emg')
 
 linkaxes(ax2, 'x');
  xlim([T(1) T(end)])

 %% 
 
  figure(1)
 
 ax(1) =subplot(9,1,1:3);
 P_p = 10*log10(P);
 imagesc(T,F,P_p)
 axis xy
 
 ax(2) =subplot(9,1,4);
 plot(T, emg_binned);
 
 
 ax(3) =subplot(9,1,5);
 plot(T, delta_binned);
 
 
 ax(4) =subplot(9,1,6);
 plot(T, delta_binned./wide_binned);
 
 ax(5) =subplot(9,1,7);
 plot(T, theta_binned);
 
 ax(6) =subplot(9,1,8);
 plot(T, theta_binned./(delta_binned./wide_binned));
 
 ax(7) =subplot(9,1,9);
 plot(T, theta_binned./(delta_binned));

 linkaxes(ax, 'x')
 xlim([T(1) T(end)])

 %% short fourier transform
 f_d = [1 4];
 f_t = [6 9];
 f_w = [1 50]; 
 
%  [S, F, T] = stft(csc.data(2,:), csc.cfg.hdr{1}.SamplingFrequency, 'Window', hamming(csc.cfg.hdr{1}.SamplingFrequency*40));
 
 [~, F, T,P] =  spectrogram(csc.data(1,:),rectwin(csc.cfg.hdr{1}.SamplingFrequency*8),csc.cfg.hdr{1}.SamplingFrequency*4,0.5:.5:50,csc.cfg.hdr{1}.SamplingFrequency);
 %%
 P_p = 10*log10(P); 
 P_p = P_p.^2;
 P_p = P_p./max(max(P_p)); 
 figure
 ax(1) = subplot(9,1,1:4);
 imagesc(T,F,P./max(max(P)))
 axis xy
 
         for iS = size(P_p,2):-1:1
           d_p(iS) = mean(P_p(find(F == f_d(1)):find(F == f_d(2)),iS));
           w_p(iS) = mean(P_p(find(F == f_w(1)):find(F == f_w(2)),iS));
           t_p(iS) = mean(P_p(find(F == f_t(1)):find(F == f_t(2)),iS));
         end
         d_norm = d_p./w_p; 
         t_norm = t_p./w_p; 
         
         td_norm = t_norm./d_norm;
 
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
 plot(T, t_p); 
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
  
 
 
 
 