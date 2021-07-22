function [pREM_idx, pREM_times, pREM_IV] =  MS_get_pREM(raw_csc, idx_in, min_len, emg_in, plot_flag)
%% MS_get_pREM: isolates phasic REM events using the methods outlined in Mizuseki et al. 2011.
%
%
%
%    Inputs:
%    - raw_csc: [struct] the 'csc struct' from MS_load_CSC (vandermeerlab
%    codebase + CEH2 function) which contains the tvec (time vector) as
%    well as data vectors for each channel and NLX headers.  this format is
%    required for filtering using FilterLFP (vandermeerlab function).
%
%    - idx_in [1 x nSamples] binary array containing the REM samples (same
%    length as raw_csc.tvec/data.
%
%    min_len [double]  value for the minimum length of a candidate event.
%    Default is 0.9sec from Mizuseki
%
%    - emg_in [1 x nSamples] an EMG signal for plotting. power from abs of
%    30-300Hz filtered EMG works well.
%
%    - plot_flag [binary] set to 1 to plot checks, set to 0 to skip.
%
%    Outputs:
%    - pREM_idx: [2 x nREM events] start and stop indices for the pREM
%    blocks.
%
%    - pREM_times: [2 x nREM events] start and stop times
%
%    - pREM_IV   [struct]  IV (interval data) which can make use of
%    vandermeerlab functions.
%
% EC 2021-06-29   initial version
%
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
%
%% initialize

if nargin < 3
    min_len = 0.9; % min duration of event
    emg_in = []; %emg for plotting.
    plot_flag = 0; % disable plots unless specified.
elseif nargin <4
    emg_in = []; %emg for plotting.
elseif nargin < 5
    plot_flag = 0; % disable plots unless specified.
end

c_ord = linspecer(5);
%% filter raw data with the Mizuseki 2011 config.
Fs = raw_csc.cfg.hdr{1}.SamplingFrequency;

% theta filter
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 11]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool

theta_csc = FilterLFP(cfg_filt_t, raw_csc); % filter the raw LFP using

theta_amp = abs(hilbert(theta_csc.data)); % get the amplitude
theta_phi = angle(hilbert(theta_csc.data)); % get the phase.
%% get REM periods

if plot_flag
    figure
    subplot(2,1,1)
end
[REM_evts, REM_IV] = MS_get_events(idx_in, plot_flag); % get the start and stop of each REM event using the REM label idx

% convert idx in IV to times
REM_IV.tstart = raw_csc.tvec(REM_IV.tstart);
REM_IV.tend = raw_csc.tvec(REM_IV.tend);

fprintf('<strong>%s</strong>: detected %d events\n', mfilename, length(REM_evts))

%% have a look at the REM events
if plot_flag
    subplot(2,1,2)
    cfg_plot = [];
    cfg_plot.display = 'tsd';
    cfg_plot.target = 'CSC6.ncs';
    PlotTSDfromIV(cfg_plot, REM_IV, raw_csc)
    xlim([raw_csc.tvec(1) raw_csc.tvec(end)])
end
%% convert REM to episode blocks

for iB = length(REM_evts):-1:1
    REM_blocks{iB} = theta_csc.data(REM_evts(iB,1):REM_evts(iB,2));
    REM_tvecs{iB} = theta_csc.tvec(REM_evts(iB,1):REM_evts(iB,2));
    REM_amp{iB} = theta_amp(REM_evts(iB,1):REM_evts(iB,2));
    REM_phi{iB} = theta_phi(REM_evts(iB,1):REM_evts(iB,2));
    
    REM_raw{iB} = raw_csc.data(REM_evts(iB,1):REM_evts(iB,2));
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

for iB = length(IPI_vec):-1:1 % still working with blocks rather than concatenated values to avoid start/end overlap.
    
    IPI_idx = IPI_vec{iB} < L5_prctile; % crit 2 get blocks < 5th prctile of smoothed IPI
    
    % Crit 1: keep blocks that are longer than 900ms
    Phasic_blocks   = MS_get_events(IPI_idx);
    
    block_dur =(Phasic_blocks(:,2) - Phasic_blocks(:,1))/Fs; % get block duration and convert to time in S.
    
    keep_blocks = block_dur > min_len; % keep only blocks that are
    
    Phasic_blocks(~keep_blocks,:) = [];
    
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
for ii  = size(pREM_times,1):-1:1
    pREM_idx(ii,1) = find(raw_csc.tvec == pREM_times(ii,1));
    pREM_idx(ii,2) = find(raw_csc.tvec == pREM_times(ii,2));
end

if isempty(pREM_times)
    fprintf('<strong>%s:  No pREM events detected</strong>\n', mfilename)
    save(['No pREM events detected ' date '.txt']);
    pREM_idx = []; pREM_IV = [];
    
    return
end
% create an IV struct
pREM_IV = iv([pREM_times(:,1), pREM_times(:,2)]);

%% plot some segments.
if plot_flag
    win_s = 2; % add some extra data
    
    for iR =1:size(pREM_idx,1)
        
        Phasic_data{iR} = raw_csc.data((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
        
        if ~isempty(emg_in) % prep EMG if you have it.
            Phasic_EMG{iR} = emg_in((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
        end
        
        figure(iR)
        cwt(Phasic_data{iR}, Fs);
        x_lim = xlim;
        hold on
        xline(win_s, '--k', 'start', 'linewidth', 2)
        xline(x_lim(2) - win_s, '--k', 'start', 'linewidth', 2)
        yline(10, '--w', '10hz', 'linewidth', 2)
        
        title(['REM event #' num2str(iR) ])
        AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(Phasic_data{iR}),Fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
        ylim([1 140])
        
        %add the LFP
        temp_tvec =  raw_csc.tvec((pREM_idx(iR,1)- win_s*Fs):(pREM_idx(iR,2)+ win_s*Fs));
        plot(temp_tvec - temp_tvec(1), (Phasic_data{iR}*1200)+4, 'w')
        
        % add EMG if you have it.
        if ~isempty(emg_in)
            plot(temp_tvec - temp_tvec(1), (Phasic_EMG{iR}*1200)+2, 'color', [.7 .7 .7])
        end
        pause(1)
        
    end
end

%% print some basic stats

fprintf('<strong>%s</strong>: pREM events detected totalling %d seconds (%.2f %% of REM)\n',...
    mfilename, numel(Phasic_data), (sum(cellfun('length',Phasic_data))/sum(cellfun('length', REM_blocks)))*100)
