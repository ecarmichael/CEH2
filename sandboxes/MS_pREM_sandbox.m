%% phasic rem detector sandbox 2021-06-03

%Follows the mizuseki et al. 2011 method: 
    % 'Detection of phasic REM. REM epochs were detected as described above.
    % To detect phasic REM epochs, we first band-pass filtered (5–12 Hz) LFP traces
    % during REM epochs, yielding y(t). The amplitudes of theta oscillations were
    % derived from Hilbert transform of y(t), and peaks of theta oscillations were
    % detected as the positive-to-negative zero crossings of the derivative dy/dt.
    % Interpeak intervals were smoothed using an 11-sample rectangular kernel.
    % Candidate epochs were detected if smoothed interpeak intervals were shorter
    % than the 10th percentile of smoothed interpeak intervals. The candidate epochs
    % were identified as phasic REM epochs if the following criteria were all fulfilled.
    % First, the duration of an epoch was longer than 900 ms. Second, the minimum of
    % smoothed interpeak intervals during an epoch was shorter than the 5th percentile
    %
    % of smoothed interspike intervals. Third, the mean amplitude of theta oscilla-
    % tions during an epoch was larger than the mean amplitude of theta oscillations
    %
    % during the entire REM sleep. A total of 5,844 s (3.68 % of REM sleep episodes)
    % was identified as phasic REM epochs.'
    
    
% criteria 
%     1. >900ms duration
%     2. epoch has min smoothed IPI < 5th percentile of all IPI_smoothed
%     3. epoch theta amp > mean theta amp for all REM. 


LFP_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo'; 
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
cd('C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch')


c_ord = linspecer(5); % just some nicer colours. 

%% load some pre-cut data

load('ms_resize.mat');

inter_dir = cd; % where to write back the sleep scored data.  
%% get the session name
parts = strsplit(cd,  filesep);

session = parts{end};
subject = parts{end-1};
date = parts{end}(1:10); 
type = strsplit(parts{end}, '_');
type = type{end}; 
if strcmp(type,'HATSwitch') && contains(subject, 'PV1069')
    type = 'HATD6_switch'; 
end

%%  manually score some data

% get the LFP files
cd(LFP_dir)

this_LFP_dir = MS_list_dir_names(cd, {subject, type}); 

cd(this_LFP_dir{1});

% get the channels to load from the pre-cut data. First channel should be
% EMG and second is the best LFP. 
cfg_load = [];
for iC = 1:length(ms_seg_resize.NLX_csc{1}.cfg.hdr) % loop over channels. 
    cfg_load.fc{iC} = [ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.AcqEntName '.ncs'];
    cfg_load.desired_sampling_frequency = ms_seg_resize.NLX_csc{1}.cfg.hdr{iC}.SamplingFrequency; 
end

% load some data. 
CSC = MS_LoadCSC(cfg_load); 


%% score the sleep data
 emg_h = abs(hilbert(CSC.data(1,:))); % get the emg power for plotting.        
 
 cfg_sleep = [];
 cfg_sleep.tvec_range = [0 5];  % number of seconds per window.
 cfg_sleep.emg_range = [min(emg_h) mean(emg_h) + std(emg_h)*5]; % default, should be based on real data.
 cfg_sleep.emg_chan = 1; % emg channel.  Can be empty.
 cfg_sleep.lfp_chans = 1; % lfp channels to be plotted can be empty. Can be 1 or more, but best to keep it less than 3. should be rows in csc.data.
 cfg_sleep.state_name = {'Wake',       'SWS',       'REM',    'Quiescence','Transition','pREM',  'Redo',     'Exit'}; %
 cfg_sleep.state_keys = {'rightarrow','uparrow', 'downarrow', 'leftarrow', 'numpad0',   'numpad1' 'backspace','backquote' }; % which key to press for each state_val
 

 score = MS_Sleep_score_UI(cfg_sleep, CSC.tvec,CSC.data(2,:), emg_h);



%% write the hypno back to the intermediate dir.  



%% grab the REM blocks and filter them with the Mizuseki 2011 config. 

REM_idx = find(contains(ms_seg_resize.hypno_label, 'REM'));

cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 12]; % to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  .

for iR = REM_idx
    % extract raw LFP only.
    this_csc = ms_seg_resize.NLX_csc{iR};
    this_csc.data = this_csc.data(2,:);
    this_csc.label = this_csc.label{2};
    temp_hdr = this_csc.cfg.hdr{2};
    this_csc.cfg.hdr = [];
    this_csc.cfg.hdr{1} = temp_hdr;
    
    theta_csc{iR} = FilterLFP(cfg_filt_t, this_csc);
    
end

%% extract the inter peak interval 

for iR = REM_idx % loop over rem episodes. 


% get the amplitude

amp{iR} = abs(hilbert(theta_csc{iR}.data)); 

Phi{iR} = angle(hilbert(theta_csc{iR}.data)); 

% get the negative to positive crossings in the phase (peaks and troughs
% are all relative)
peaks = []; 
for ii = 1:length(Phi{iR})-1
   if (Phi{iR}(ii)>0) && (Phi{iR}(ii+1)<=0)
    peaks = [peaks ii+1]; 
   end
end

% get the Inter Peak Interval
IPI{iR} = diff(theta_csc{iR}.tvec(peaks));

% smooth with 11 sample rec window.  warning, gives much lower IPIs. Use
% for distribution only. 
IPI_smooth{iR} = conv2(IPI{iR},rectwin(11), 'same'); 



% [TODO] fill in IPI values per cycle to match the actual data. 






% make a sample plot if needed
% figure(101)
% ax(1) = subplot(3,1,1);
% hold on
% plot(theta_csc{iR}.tvec, ms_seg_resize.NLX_csc{iR}.data(2,:), 'k')
% plot(theta_csc{iR}.tvec, theta_csc{iR}.data, 'b')
% plot(theta_csc{iR}.tvec, amp{iR}, 'g')
% plot(theta_csc{iR}.tvec(peaks), theta_csc{iR}.data(peaks), 'x', 'color', 'r')
% 
% 
% ax(2) = subplot(3,1,2);
% hold on
% plot(theta_csc{iR}.tvec, Phi{iR}, 'r')
% plot(theta_csc{iR}.tvec(2:end), diff(Phi{iR}), '--r')
% plot(theta_csc{iR}.tvec(peaks), Phi(peaks), 'x', 'color', 'r')
% 
% linkaxes(ax, 'x')
% 
% subplot(3,1,3)
% hold on
% histogram(IPI,25)
% histogram(IPI_smooth,25)
% legend('IPI', 'IPI smoothed')

% figure(102)
% ax(1) = subplot(2,2,1:2);
% hold on
% plot(theta_csc{iR}.tvec, theta_csc{iR}.data,'color',  c_ord(1,:))
% plot(theta_csc{iR}.tvec, amp{iR}, 'color', c_ord(2,:))
% ax(2) = subplot(2,2,3:4);
% hold on
% plot(theta_csc{iR}.tvec(peaks(1:end-1)), IPI_smooth{iR}, 'color',  c_ord(4,:))
% plot(theta_csc{iR}.tvec(peaks(1:end-1)), IPI{iR}, 'color', c_ord(1,:))
% linkaxes(ax, 'x')
% 
% pause
% close all

end
%%  collect the IPIs to get a distribution
all_IPI = []; all_IPI_smooth = []; % collect the ISI for crit 1
all_amp = []; % collect the amplitude for crit 3

for iR = REM_idx
    all_IPI = [all_IPI; IPI{iR}]; 
    
    all_IPI_smooth = [all_IPI_smooth; IPI_smooth{iR}]; 
    
    all_amp = [all_amp, amp{iR}];
end

L10_prctile = prctile(all_IPI_smooth, 10); 
L5_prctile = prctile(all_IPI_smooth, 5); 
L50_prctile = prctile(all_IPI_smooth, 50);


subplot(3,2,5)
histogram(all_IPI,50, 'facecolor', c_ord(1,:));
xlabel('IPI')
legend('IPI')

subplot(3,2,6)
histogram(all_IPI_smooth,50,'facecolor', c_ord(4,:));
vline([L10_prctile,L5_prctile, L50_prctile], {'k', 'r', 'm'}, {'L10', 'L5', 'L50'}); 
legend('IPI smoothed')
xlabel('IPI')


%% Crit 1 find blocks with >900ms duration. 



%% Crit 2 remove blocks without a min IPI smoothed < 5th percentile of all IPI smoothed


%% crit 3 Remove blocks with mean theta amp < mean theta amplitude for all REM. 


mean_theta_amp = mean(all_amp); 


