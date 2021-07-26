%% phasic rem detector sandbox 2021-06-03

%Follows the mizuseki et al. 2011 method:
% 'Detection of phasic REM. REM epochs were detected as described above.
% To detect phasic REM epochs, we first band-pass filtered (5�12 Hz) LFP traces
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

if isunix
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'));
    addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
    data_dir = '/home/williamslab/Desktop/Jisoo_sleep_LFP/PV1060/11_26_2019_PV1060_HATSwitch';
    LFP_dir = '/home/williamslab/Desktop/Jisoo_sleep_LFP';
else
    
    LFP_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo';
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
   data_dir ='J:\Williams_Lab\Jisoo\Jisoo_Project\Inter\PV1069\10_18_2019_PV1069_HATD5'; % change this to the data folder that you want.  LFP will update automatically.
   cell_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\4.PlaceCell'; % where to find place cell classification and centroids. 
   
   % data_dir = 'C:\Users\ecarm\Dropbox (Williams
   % Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch';
   % my drop box doesn't work here and whole data is in J
end

cd(data_dir)
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

%%  prepare LFP data

% get the LFP files
cd(LFP_dir)

this_LFP_dir = MS_list_dir_names(cd, {subject, type});
% this_LFP_dir{1} = '/home/williamslab/Desktop/Jisoo_sleep_LFP/2019-07-15_09-50-03_PV1060_LTD1';

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
EVT = LoadEvents([]); % load the events file.


%% plot a quick PSD check

% MS_Quick_psd()

%% prepare the data for scoring.
% since JC data has a gap in a single recording for the track, account for
% that here.
S_rec_idx = find(contains(EVT.label, 'Starting Recording')); % get the index for the start recoding cell
if length(EVT.t{S_rec_idx}) >1
    fprintf('Multiple (%d) Recordings in one continuous CSC.  appending for sleep scoring spec\n', numel(EVT.t{S_rec_idx}))
    
    all_tvec = [];
    for iCut = 2:numel(EVT.t{S_rec_idx})
        if iCut == numel(EVT.t{S_rec_idx})
            tvec_cut = CSC.tvec(nearest_idx(EVT.t{S_rec_idx}(iCut), CSC.tvec):end);
        else
            tvec_cut = CSC.tvec(nearest_idx(EVT.t{S_rec_idx}(iCut), CSC.tvec):iCut+1);
        end
        
        tvec_cut = tvec_cut - tvec_cut(1);
        all_tvec = [all_tvec, tvec_cut];
    end
    start_tvec = CSC.tvec(1:nearest_idx(EVT.t{S_rec_idx+1}(1),CSC.tvec));
    all_tvec = [start_tvec; all_tvec+start_tvec(end)+(1/CSC.cfg.hdr{1}.SamplingFrequency)];
    
    % update CSC
    CSC_cut = CSC;
    CSC_cut.tvec = all_tvec;
end

% extract the EMG
CSC_emg = CSC_cut;
CSC_emg.data = CSC_emg.data(1,:);
CSC_emg.cfg.hdr = [];
CSC_emg.cfg.hdr{1} = CSC_cut.cfg.hdr{1};


% filter the EMG
cfg_emg = [];
cfg_emg.f = [15 300];
cfg_emg.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_emg.order = 16; %type filter order
emg_f = FilterLFP(cfg_emg, CSC_emg);

emg_h = abs(hilbert(emg_f.data)); % get the emg power for plotting.

emg_rms = sqrt(movmean(emg_f.data.^2, 1000));                         % RMS Value Over ‘WinLen’ Samples

%% score the sleep data

cfg_sleep = [];
cfg_sleep.tvec_range = [0 5];  % number of seconds per window.
cfg_sleep.emg_range = [min(emg_h) mean(emg_h) + std(emg_h)*5]; % default, should be based on real data.
cfg_sleep.emg_chan = 1; % emg channel.  Can be empty.
cfg_sleep.lfp_chans = 1; % lfp channels to be plotted can be empty. Can be 1 or more, but best to keep it less than 3. should be rows in csc.data.
cfg_sleep.state_name = {'Wake',       'SWS',       'REM',    'Quiescence','Transition','pREM',  'Redo',     'Exit'}; %
cfg_sleep.state_keys = {'rightarrow','uparrow', 'downarrow', 'leftarrow', 't',   'numpad1' 'backspace','backquote' }; % which key to press for each state_val


score = MS_Sleep_score_UI(cfg_sleep, CSC_cut.tvec,CSC.data(2,:), emg_h);



%% write the hypno back to the intermediate dir.

cd(inter_dir); % commented by jisoo maybe add creating folder here?
mkdir('pREM')

if exist('score', 'var')
    % put it in a nice format.
    Hypno = [];
    Hypno.tvec = CSC_cut.tvec;
    Hypno.data = score;
    Hypno.labels = cfg_sleep.state_name;
    Hypno.cfg.date = date;
    Hypno.cfg.method = 'MS_Sleep_score_UI';
    Hypno.cfg.notes = 'pREM was not scored in mannual screening';
    
    % save it.
    cd(inter_dir);
    save([inter_dir filesep 'pREM' filesep 'Hypno.mat'], 'Hypno', '-v7.3')
    save('Hypno.mat', 'Hypno', '-v7.3')
    
else
    fprintf('<strong>No Score var found. skipping saving</strong>\n')
end
%% filter raw data with the Mizuseki 2011 config.


cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 12]; % to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.

% extract raw LFP only.
this_csc =CSC_cut;
this_csc.data = this_csc.data(2,:);
this_csc.label = this_csc.label{2};
temp_hdr = this_csc.cfg.hdr{2};
this_csc.cfg.hdr = [];
this_csc.cfg.hdr{1} = temp_hdr;

theta_csc = FilterLFP(cfg_filt_t, this_csc);

theta_amp = abs(hilbert(theta_csc.data));
theta_phi = angle(hilbert(theta_csc.data));
%% get REM periods

if ~exist('Hypno', 'var')
    if exist([inter_dir filesep 'Hypno.mat'], 'file')
        fprintf('<strong>Loading Hypno from </strong> %s\n', inter_dir);
        load([inter_dir filesep 'Hypno.mat']);
    end
else
    fprintf('Hypno present. Using var\n');
end

REM_label_idx = find(contains(Hypno.labels, 'REM'));
REM_idx = Hypno.data == REM_label_idx;
REM_idx = sum(REM_idx,2); % if there are more than one REM label (ie REM & pREM) keep both.
REM_idx(REM_idx >0) = 1;  % in case there is any overlap.
REM_idx = logical(REM_idx); % make it a logical again.

figure
subplot(2,1,1)
[REM_evts, REM_IV] = MS_get_events(REM_idx', 1); % get the start and stop of each REM event using the REM label idx

% convert idx in IV to times
REM_IV.tstart = CSC_cut.tvec(REM_IV.tstart);
REM_IV.tend = CSC_cut.tvec(REM_IV.tend);


%% have a look at the REM events
subplot(2,1,2)
cfg_plot = [];
cfg_plot.display = 'tsd';
cfg_plot.target = 'CSC6.ncs';
PlotTSDfromIV(cfg_plot, REM_IV, CSC_cut)
xlim([CSC_cut.tvec(1) CSC_cut.tvec(end)])

%% convert REM to episode blocks

for iB = length(REM_evts):-1:1
    REM_blocks{iB} = theta_csc.data(REM_evts(iB,1):REM_evts(iB,2));
    REM_tvecs{iB} = theta_csc.tvec(REM_evts(iB,1):REM_evts(iB,2));
    REM_amp{iB} = theta_amp(REM_evts(iB,1):REM_evts(iB,2));
    REM_phi{iB} = theta_phi(REM_evts(iB,1):REM_evts(iB,2));
    
    REM_raw{iB} = this_csc.data(REM_evts(iB,1):REM_evts(iB,2));
end

%% extract the inter peak interval

for iB = length(REM_blocks):-1:1
    % get the negative to positive crossings in the phase (peaks and troughs
    % are all relative)
    phi_peaks = [];
    for ii = 1:length(REM_phi{iB})-1
        if (REM_phi{iB}(ii)>0) && (REM_phi{iB}(ii+1)<=0)
            phi_peaks = [phi_peaks ii+1];
        end
    end
    
    % get the Inter Peak Interval
    IPI{iB} = diff(REM_tvecs{iB}(phi_peaks));
    
    % smooth with 11 sample rec window.  warning, gives much lower IPIs. Use
    % for distribution only.
    IPI_smooth{iB} = conv2(IPI{iB},rectwin(11), 'same');
    
    % [TODO] fill in IPI values per cycle to match the actual data.
    IPI_vec{iB} = nan(size(REM_phi{iB}));
    
    for ii = 1:length(phi_peaks)
        if ii == 1 % fill in the first data point to the first peak
            IPI_vec{iB}(1:phi_peaks(ii+1)) = IPI_smooth{iB}(ii);
        elseif ii == length(phi_peaks)
            IPI_vec{iB}(phi_peaks(ii):end) = IPI_smooth{iB}(ii-1);
        else
            IPI_vec{iB}(phi_peaks(ii):phi_peaks(ii+1)) = IPI_smooth{iB}(ii);
        end
    end
    
    
    
    % make a sample plot if needed
    figure(iB)
    ax(1) = subplot(3,2,1:2);
    hold on
    yyaxis left
    plot(REM_tvecs{iB}, REM_blocks{iB}, 'k')
    plot(REM_tvecs{iB}, REM_amp{iB},'--', 'color', c_ord(3,:), 'linewidth', 2)
    plot(REM_tvecs{iB}(phi_peaks), REM_blocks{iB}(phi_peaks), 'x', 'color', c_ord(2,:))
    ylabel('voltage');
    
    yyaxis right
    plot(REM_tvecs{iB}, IPI_vec{iB}, 'color', c_ord(5,:), 'linewidth', 2);
    ylabel('Smoothed IPI');
    
    ayy = gca;
    ayy.YAxis(1).Color = 'k';
    ayy.YAxis(2).Color = c_ord(5,:);
    
    ax(2) = subplot(3,2,3:4);
    hold on
    plot(REM_tvecs{iB}, REM_phi{iB}, 'k')
    plot(REM_tvecs{iB}(2:end), diff(REM_phi{iB}), '--r')
    plot(REM_tvecs{iB}(phi_peaks), REM_phi{iB}(phi_peaks), 'x', 'color', 'r')
    
    linkaxes(ax, 'x')
    
    subplot(3,2,5)
    histogram(IPI{iB},25,'facecolor', c_ord(1,:));
    legend('IPI')
    x_label = get(gca, 'xtick');
    % set(gca, 'xticklabel', round((1./x_label)*100)/100)
    for ii = 1:length(x_label)
        vline(x_label(ii), '--k', {num2str(round((1./x_label(ii))*100)/100)});
    end
    xlabel('IPI (s)')
    
    subplot(3,2,6)
    histogram(IPI_smooth{iB},25,'facecolor', c_ord(4,:));
    legend('IPI smoothed (s)')
    
    % figure(102)
    % ax(1) = subplot(2,2,1:2);
    % hold on
    % plot(theta_csc{iR}.tvec, theta_csc{iR}.data,'color',  c_ord(1,:))
    % plot(theta_csc{iR}.tvec, amp{iR}, 'color', c_ord(2,:))
    % ax(2) = subplot(2,2,3:4);
    % hold on
    % plot(theta_csc{iR}.tvec(phi_peaks(1:end-1)), IPI_smooth{iR}, 'color',  c_ord(4,:))
    % plot(theta_csc{iR}.tvec(phi_peaks(1:end-1)), IPI{iR}, 'color', c_ord(1,:))
    % linkaxes(ax, 'x')
    
    pause(.25)
    % close all
end


close all
%%  collect the IPIs to get a distribution
all_IPI = []; all_IPI_smooth = []; % collect the ISI for crit 1
all_amp = []; % collect the amplitude for crit 3

for iB = 1:length(REM_blocks)
    all_IPI = [all_IPI; IPI{iB}];
    
    all_IPI_smooth = [all_IPI_smooth; IPI_smooth{iB}];
    
    all_amp = [all_amp, REM_amp{iB}];
end

L10_prctile = prctile(all_IPI_smooth, 10);
L5_prctile = prctile(all_IPI_smooth, 5);
L50_prctile = prctile(all_IPI_smooth, 50);

figure(101)
subplot(1,2,1)
histogram(all_IPI,50, 'facecolor', c_ord(1,:));
x_label = get(gca, 'xtick');
% set(gca, 'xticklabel', round((1./x_label)*100)/100)
for ii = 1:length(x_label)
    vline(x_label(ii), '--k', {num2str(round((1./x_label(ii))*100)/100)});
end
xlabel('IPI (s)')
xlabel('IPI')
legend('IPI')

subplot(1,2,2)
histogram(all_IPI_smooth,50,'facecolor', c_ord(4,:));
vline([L10_prctile,L5_prctile, L50_prctile], {'k', 'r', 'm'}, {'L10', 'L5', 'L50'});
legend('IPI smoothed')
xlabel('IPI')


%% Crit 2 remove blocks without a min IPI smoothed < 5th percentile of all IPI smoothed

min_len = .9;

for iB = length(IPI_vec):-1:1 % still working with blocks rather than concatenated values to avoid start/end overlap.
    
    IPI_idx = IPI_vec{iB} < L5_prctile; % crit 2 get blocks < 5th prctile of smoothed IPI
    
    % Crit 1: keep blocks that are longer than 900ms
    Phasic_blocks   = MS_get_events(IPI_idx);
    
    block_dur =(Phasic_blocks(:,2) - Phasic_blocks(:,1))/CSC_cut.cfg.hdr{1}.SamplingFrequency; % get block duration and convert to time in S.
    
    keep_blocks = block_dur > min_len; % keep only blocks that are
    
    Phasic_blocks(~keep_blocks,:) = [];
    keep_t_amp = []; 
    if ~isempty(Phasic_blocks)
        % Crit 3: theta amp in block must be great than mean of all theta
        for ii = 1:size(Phasic_blocks,1)
            if mean(REM_amp{iB}(Phasic_blocks(ii,1):Phasic_blocks(ii,2)))  > mean(theta_amp)
                keep_t_amp(ii) = 1;
            else
                keep_t_amp(ii) = 0;
            end
        end
        
        Phasic_blocks(~keep_t_amp,:) = []; % remove low theta blocks.
    end
    
    %convert back to the time vector.
    Phasic_times{iB} = REM_tvecs{iB}(Phasic_blocks);
    
    Phasic_idx{iB} = Phasic_blocks;
    
end

%% compile events and check with some plots
pREM_times = []; pREM_idx = [];

for iB = 1:length(Phasic_times)
    if isempty(Phasic_times{iB})
        continue
    end
    if size(Phasic_times{iB},2) > 1 % hack due to some dimension flipping taking place.  There is a simple one line solution that I am missing.
        pREM_times = [pREM_times; Phasic_times{iB}];
    else
        pREM_times = [pREM_times; Phasic_times{iB}'];
    end
end


% convert to indexes
for ii  = 1:size(pREM_times,1)
    pREM_idx(ii,1) = find(this_csc.tvec == pREM_times(ii,1));
    pREM_idx(ii,2) = find(this_csc.tvec == pREM_times(ii,2));
end



%% plot some segments.
win_s = 2; % add some extra data
Fs = this_csc.cfg.hdr{1}.SamplingFrequency; % sampling freq

% save plots.
cd(inter_dir);
mkdir('pREM')

for iR =1:size(pREM_idx,1)
    
    Phasic_data{iR} = this_csc.data((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
    Phasic_EMG{iR} = emg_h((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
    
    cwt(Phasic_data{iR}, Fs);
    x_lim = xlim;
    hold on
    xline(win_s, '--k', 'start', 'linewidth', 2);
    xline(x_lim(2) - win_s, '--k', 'start', 'linewidth', 2);
    yline(10, '--w', '10hz', 'linewidth', 2);
    
    title(['REM event #' num2str(iR) ]);
    AX = gca;
    [minf,maxf] = cwtfreqbounds(numel(Phasic_data{iR}),Fs);
    
    freq = 2.^(round(log2(minf)):round(log2(maxf)));
    AX.YTickLabelMode = 'auto';
    AX.YTick = freq;
    ylim([1 140]);
    
    %add the LFP
    temp_tvec =  this_csc.tvec((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
    plot(temp_tvec - temp_tvec(1), (Phasic_data{iR}*1200)+4, 'w');
    
    plot(temp_tvec - temp_tvec(1), (Phasic_EMG{iR}*1200)+2, 'color', [.7 .7 .7]);
    
    %         pause(1)
    
    saveas(gcf,[inter_dir filesep 'pREM' filesep 'pREM_event_' num2str(iR) '.png'])
end


%% export the intervals and times

save([inter_dir filesep 'pREM' filesep 'pREM_idx.mat'], 'pREM_idx');
save([inter_dir filesep 'pREM' filesep 'pREM_times.mat'], 'pREM_times');

fprintf('<strong>%s</strong>: pREM events detected totalling %d seconds (%.2f %% of REM)\n',...
    mfilename, numel(Phasic_data), (sum(cellfun('length',Phasic_data))/sum(cellfun('length', REM_blocks)))*100)


%% get the Ca time and convert the pREM episodes to that time. 
load('ms_resize.mat')


%% align the CA and the NLX times. 
all_aligned_tvec_pre = []; 
all_aligned_tvec_post = []; 
all_evts = []; 

for iSeg = 1:length(ms_seg_resize.time)
%   if length(ms_seg_resize.time{iSeg}) ~= length(ms_seg_resize.NLX_evt{iSeg}.t{end})
%         fprintf('Segment # %s, length of time %d and NLX_events  %d do not match. Using ms_seg.time to generate NLX TS\n',  num2str(iSeg),length(ms_seg_resize.time{iSeg}),length(ms_seg_resize.NLX_evt{iSeg}.t{end}))
%         c_ms_time = ((ms_seg_resize.time{iSeg}-ms_seg_resize.time{iSeg}(1))*0.001)';
%         %         c_nlx_time = ((ms_seg.NLX_evt.t{end}-ms_seg.NLX_evt.t{end}(1)))%-0.0077;
%         
%         %         NLX_ts = [c_nlx_time, c_ms_time(length(c_nlx_time)+1:end)]+ms_seg.NLX_evt.t{end}(1);
%         
%         NLX_ts = c_ms_time + ms_seg_resize.NLX_evt{iSeg}.t{end}(1);
%     else
        NLX_ts = ms_seg_resize.NLX_evt{iSeg}.t{end};
        
%   end
  
  all_evts = [all_evts NLX_ts]; 
  
  if strcmpi(ms_seg_resize.pre_post{iSeg}, 'pre')
    all_aligned_tvec_pre = [all_aligned_tvec_pre, NLX_ts];
  elseif strcmpi(ms_seg_resize.pre_post{iSeg}, 'post')
    all_aligned_tvec_post = [all_aligned_tvec_post, NLX_ts];
  end
end

% get the index for the concatenated Ca data which is split into pre and
% post blocks (from MS_re_binarize_JC.m or MS_extract_AMP_Phi.m)

all_pREM_Ca_idx = NaN(length(pREM_times),2); 
% pREM_Ca_idx_post = NaN(length(pREM_times),2);
for ii = 1:length(pREM_times)
    if (pREM_times(ii,1) < all_aligned_tvec_post(1)) || (pREM_times(ii,2) < all_aligned_tvec_post(1))
        
        all_pREM_Ca_idx(ii,1) = nearest_idx3(pREM_times(ii,1), all_aligned_tvec_pre)';
        all_pREM_Ca_idx(ii,2) = nearest_idx3(pREM_times(ii,2), all_aligned_tvec_pre)';
        
        % catch for pREM that fall outside of the Ca recording blocks.
        if all_pREM_Ca_idx(ii,2) - all_pREM_Ca_idx(ii,1) <=0
            all_pREM_Ca_idx(ii,:) = NaN;
        end
        pREM_block{ii} = 'pre';

    else
        all_pREM_Ca_idx(ii,1) = nearest_idx3(pREM_times(ii,1), all_aligned_tvec_post)';
        all_pREM_Ca_idx(ii,2) = nearest_idx3(pREM_times(ii,2), all_aligned_tvec_post)';
        
        % catch for pREM that fall outside of the Ca recording blocks.
        if all_pREM_Ca_idx(ii,2) - all_pREM_Ca_idx(ii,1) <=0
            all_pREM_Ca_idx(ii,:) = NaN ;
        end
                pREM_block{ii} = 'post';
    end

end

all_pREM_CA.idx = all_pREM_Ca_idx; 
all_pREM_CA.label = pREM_block; 


%% save and print number of pREM events that overlap with Ca recording. 
save([inter_dir filesep 'pREM' filesep 'pREM_CA_idx.mat'], 'all_pREM_CA');

fprintf('<strong>%s</strong>:  %d/%d pREM events occured during Miniscope recording. \n',...
    mfilename, sum(~isnan(all_pREM_Ca_idx(:,1))), size(all_pREM_Ca_idx,1))

% all_pREM_Ca_idx = [all_pREM_Ca_idx; pREM_Ca_idx_post];

% % remove zero values (space holder for pREM's outside of Ca recordings. 
% pREM_Ca_idx_pre = pREM_Ca_idx_pre(any(pREM_Ca_idx_pre,2),:); 
% pREM_Ca_idx_post = pREM_Ca_idx_post(any(pREM_Ca_idx_post,2),:); 
% 

% % test plot
% figure(101)
% hold on
% plot(diff(all_aligned_tvec_post));
%  plot(pREM_Ca_idx_post(:,1)-1,50*ones(size(pREM_Ca_idx_post(:,1))), 'xr')
%   plot(pREM_Ca_idx_post(:,2)-1,50*ones(size(pREM_Ca_idx_post(:,2))), 'xg')

%% 
figure(1010)
hold on
plot(CSC_cut.tvec, CSC_cut.data(2,:))
for ii = 1:length(pREM_times)
    xline(CSC_cut.tvec(pREM_idx(ii,1)), 'b');
    xline(CSC_cut.tvec(pREM_idx(ii, 2)), 'r');
end
plot(all_evts, ones(size(all_evts))*median(CSC_cut.data(2,:)), 'o')
vline(EVT.t{1}(2), '--r', '--> post')

%% make some plots with corresponding rasters

% get the binarized pre and post sleep data 
load('all_binary_pre.mat');
load('all_binary_post.mat'); 



% get the cell centroids for coloring.  (if they exist)
if exist([cell_dir filesep lower(subject) filesep type filesep 'SA.mat'], 'file')
    load([cell_dir filesep lower(subject) filesep type filesep 'spatial_analysis.mat']);
    
    place_idx = zeros(length(spatial_analysis.raw),1); % allocate the index array
    centroids = nan(size(place_idx));
    
    for iC = length(spatial_analysis.raw):-1:1
        if spatial_analysis.raw{iC,3}.IsPlaceCell
            place_idx(iC) = 1;
            centroids(iC) = spatial_analysis.raw{iC,3}.PlaceFieldCentroid{1}(1);
        end
    end
    
   [~, cent_sort] = sort(centroids); 
   place_idx = place_idx(cent_sort); 

else
    cent_sort = 1:size(all_binary_post,2); % just use default sort. 
    place_idx = ones(size(cent_sort)); 
end



for ii  = 1:length(all_pREM_Ca_idx)
%     close all

    if ~isnan(all_pREM_Ca_idx(ii,1))

%         Phasic_data{iR} = this_csc.data((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
%         Phasic_EMG{iR} = emg_h((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
        cwt(Phasic_data{ii}, Fs);
        x_lim = xlim;
        hold on
        xline(win_s, '--k', 'start', 'linewidth', 2);
        xline(x_lim(2) - win_s, '--k', 'end', 'linewidth', 2);
        yline(10, '--w', '10hz', 'linewidth', 2);
        
        title(['REM event #' num2str(ii) ]);
        AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(Phasic_data{ii}),Fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
        ylim([1 140]);
        
        %add the LFP
        temp_tvec =  this_csc.tvec((pREM_idx(ii,1)- win_s*Fs):(pREM_idx(ii,2)+ win_s*Fs));
        plot(temp_tvec - temp_tvec(1), (Phasic_data{ii}*1200)+4, 'w');
        
        plot(temp_tvec - temp_tvec(1), (Phasic_EMG{ii}*1200)+2, 'color', [.7 .7 .7]);
        
        
        % make a raster
        if strcmpi(all_pREM_CA.label(ii), 'pre')
            this_ca = all_binary_pre(all_pREM_Ca_idx(ii,1) - win_s*30:all_pREM_Ca_idx(ii,2) + win_s*30,:)';
        elseif strcmpi(all_pREM_CA.label(ii), 'post')
            this_ca = all_binary_post(all_pREM_Ca_idx(ii,1) - win_s*30:all_pREM_Ca_idx(ii,2) + win_s*30,:)';
        end
        this_tvec = (all_pREM_Ca_idx(ii,1) - win_s*30:all_pREM_Ca_idx(ii,2) + win_s*30);
        this_tvec = (this_tvec - this_tvec(1))/30;
        
        figure(2);
        c_mat = [linspecer(sum(place_idx));  repmat([1 1 1], sum(~place_idx),1)]; % make colors depending on the 
        MS_Ca_Raster(this_ca(cent_sort,:),this_tvec, 6, c_mat);%repmat([1,1,1], size(this_ca,1),1)
        xline(win_s, '--w', 'start', 'linewidth', 2);
        x_lim = xlim;
        xline(x_lim(2) - win_s, '--w', 'end', 'linewidth', 2);
        set(gca, 'color', 'k'); %set background color. 
        colormap([linspecer(sum(place_idx));  repmat([1 1 1], 1,1)]); 
        cx = colorbar; 
        cx.TickLabels = cx.Ticks * max(centroids); 
        cx.Label.String = 'place cell centroid'; 
        
        % put the plots together
%         axcp = copyobj(ax, fig2);
%         set(axcp,'Position',get(ax1,'position'));
%         delete(ax1);
        figlist=get(groot,'Children');
        
        newfig=figure;
        tcl=tiledlayout(2,1);
        
        for jj= 1:numel(figlist)
            figure(figlist(jj));
            ax=gca;
            ax.Parent=tcl;
            ax.Layout.Tile=jj;
        end
        close(1);
        close(2);
        set(gcf, 'position', [662 96 758 892])
        set(gcf, 'InvertHardcopy', 'off')

        %         pause(1)
        
        saveas(gcf,[inter_dir filesep 'pREM' filesep 'pREM_event_' num2str(ii) '_' all_pREM_CA.label{ii} '_Raster.png'])
    end
    close all
end
%% Get the Sleep state stats

fprintf('<strong>%s</strong>: pREM events detected totalling %d seconds (%.2f %% of REM)\n',...
    mfilename, numel(Phasic_data), (sum(cellfun('length',Phasic_data))/sum(cellfun('length', REM_blocks)))*100)

pre_idx = contains(pREM_block, 'pre'); 
post_idx = ~pre_idx; 


all_REM_tvec = cell2mat(REM_tvecs'); 
post_tvec_idx = nearest_idx(EVT.t{1}(2),all_REM_tvec); 

pREM_dur.all = (sum(cellfun('length',Phasic_data))/length(all_REM_tvec))*100; 
pREM_dur.pre = (sum(cellfun('length',Phasic_data(pre_idx)))/length(all_REM_tvec(1:post_tvec_idx-1)))*100;
pREM_dur.post = (sum(cellfun('length',Phasic_data(post_idx)))/length(all_REM_tvec(post_tvec_idx:end)))*100;

% get mean duration of the events

%% Crit 1 find blocks with >900ms duration.
% Fs = CSC_cut.cfg.hdr{1}.SamplingFrequency;
% dur_keep_idx = ((REM_evts(:,2) - REM_evts(:,1))./Fs) > .9; % only keep blocks longer than 900ms;
%
% % recompute the all_rems without any removed blocks.
%
% all_IPI = []; all_IPI_smooth = []; % collect the ISI for crit 1
% all_amp = []; % collect the amplitude for crit 3
%
% keep_blocks = find(dur_keep_idx);
%
% for iB = 1:length(keep_blocks)
%     all_IPI = [all_IPI; IPI{keep_blocks(iB)}];
%
%     all_IPI_smooth = [all_IPI_smooth; IPI_smooth{keep_blocks(iB)}];
%
%     all_amp = [all_amp, REM_amp{keep_blocks(iB)}];
% end





%% crit 3 Remove blocks with mean theta amp < mean theta amplitude for all REM.

% mean_t = mean(all_amp);
%
% t_amp_idx = all_amp > mean_t;
%
% % make a plot to check this;
% figure(104)
% hold on
% temp_t = (0:length(all_amp)-1)./Fs;
% plot(temp_t, all_amp, 'k');
% plot(temp_t(t_amp_idx), all_amp(t_amp_idx),'.', 'color', c_ord(2,:));
%
%

% mean_theta_amp = mean(all_amp);








%% %% REM BLOCK versionextract the inter peak interval
%
% for iR = REM_idx % loop over rem episodes.
% % get the amplitude
%
% amp{iR} = abs(hilbert(theta_csc{iR}.data));
%
% Phi{iR} = angle(hilbert(theta_csc{iR}.data));
%
% % get the negative to positive crossings in the phase (peaks and troughs
% % are all relative)
% peaks = [];
% for ii = 1:length(Phi{iR})-1
%    if (Phi{iR}(ii)>0) && (Phi{iR}(ii+1)<=0)
%     peaks = [peaks ii+1];
%    end
% end
%
% % get the Inter Peak Interval
% IPI{iR} = diff(theta_csc{iR}.tvec(peaks));
%
% % smooth with 11 sample rec window.  warning, gives much lower IPIs. Use
% % for distribution only.
% IPI_smooth{iR} = conv2(IPI{iR},rectwin(11), 'same');
%
%
%
% % [TODO] fill in IPI values per cycle to match the actual data.
%
%
% % make a sample plot if needed
% % figure(101)
% % ax(1) = subplot(3,1,1);
% % hold on
% % plot(theta_csc{iR}.tvec, ms_seg_resize.NLX_csc{iR}.data(2,:), 'k')
% % plot(theta_csc{iR}.tvec, theta_csc{iR}.data, 'b')
% % plot(theta_csc{iR}.tvec, amp{iR}, 'g')
% % plot(theta_csc{iR}.tvec(peaks), theta_csc{iR}.data(peaks), 'x', 'color', 'r')
% %
% %
% % ax(2) = subplot(3,1,2);
% % hold on
% % plot(theta_csc{iR}.tvec, Phi{iR}, 'r')
% % plot(theta_csc{iR}.tvec(2:end), diff(Phi{iR}), '--r')
% % plot(theta_csc{iR}.tvec(peaks), Phi(peaks), 'x', 'color', 'r')
% %
% % linkaxes(ax, 'x')
% %
% % subplot(3,1,3)
% % hold on
% % histogram(IPI,25)
% % histogram(IPI_smooth,25)
% % legend('IPI', 'IPI smoothed')
%
% % figure(102)
% % ax(1) = subplot(2,2,1:2);
% % hold on
% % plot(theta_csc{iR}.tvec, theta_csc{iR}.data,'color',  c_ord(1,:))
% % plot(theta_csc{iR}.tvec, amp{iR}, 'color', c_ord(2,:))
% % ax(2) = subplot(2,2,3:4);
% % hold on
% % plot(theta_csc{iR}.tvec(peaks(1:end-1)), IPI_smooth{iR}, 'color',  c_ord(4,:))
% % plot(theta_csc{iR}.tvec(peaks(1:end-1)), IPI{iR}, 'color', c_ord(1,:))
% % linkaxes(ax, 'x')
% %
% % pause
% % close all
%
% end
% %%  collect the IPIs to get a distribution
% all_IPI = []; all_IPI_smooth = []; % collect the ISI for crit 1
% all_amp = []; % collect the amplitude for crit 3
%
% for iR = REM_idx
%     all_IPI = [all_IPI; IPI{iR}];
%
%     all_IPI_smooth = [all_IPI_smooth; IPI_smooth{iR}];
%
%     all_amp = [all_amp, amp{iR}];
% end
%
% L10_prctile = prctile(all_IPI_smooth, 10);
% L5_prctile = prctile(all_IPI_smooth, 5);
% L50_prctile = prctile(all_IPI_smooth, 50);
%
%
% subplot(3,2,5)
% histogram(all_IPI,50, 'facecolor', c_ord(1,:));
% xlabel('IPI')
% legend('IPI')
%
% subplot(3,2,6)
% histogram(all_IPI_smooth,50,'facecolor', c_ord(4,:));
% vline([L10_prctile,L5_prctile, L50_prctile], {'k', 'r', 'm'}, {'L10', 'L5', 'L50'});
% legend('IPI smoothed')
% xlabel('IPI')
%
%
% %% Crit 1 find blocks with >900ms duration.
%
%
%
% %% Crit 2 remove blocks without a min IPI smoothed < 5th percentile of all IPI smoothed





