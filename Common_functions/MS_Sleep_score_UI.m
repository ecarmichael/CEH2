function score = MS_Sleep_score_UI(cfg_in, tvec, data, emg, score_in)
%% MS_Sleep_score_UI: uses a figure and button presses to move through and score an LFP recording.
%
%
%
%    Inputs:
%    - cfg [struct]   configuration see the defaults below.
%
%    - tvec [1 x nSamples]  array of time points corresponding to data and
%    emg.
%
%    - data [nChan x nSamples]  contains continuous data points (LFP, Ca
%    trace, ...) corresponding to tvec
%
%    - emg [nChan x nSamples] contains corresponding emg data (or
%    something else if you like)
%
%    [OPTIONAL]
%    - score_in [1 x nSamples]   array with score data.  can be added to
%    plots for comparisons or verification of automated methods.
%
%    Outputs:
%    - score_out  [1 x nSamples]  contains score values for each timepoint.
%
%    - IV_out   [struct]  IV struct containing start and stop times for
%    each sleep category.  Great for fast plotting and restricting data.
%
%
%
%
% EC 2020-05-25   initial version
%
%
%
%% set inputs;

if nargin < 3
    error('Requires a cfg, lfp_tvec, and lfp_data')
elseif nargin == 3
    emg = [];
    score_in =[];
elseif nargin == 4
    score_in = [];
elseif nargin ==5
    score_in = [];
end


%% initialize
cfg_def = [];
cfg_def.tvec_range = [0 30];  % number of seconds per window.
cfg_def.emg_range = [-0.001 0.001]; % default, should be based on real data.
cfg_def.emg_chan = []; % emg channel.  Can be empty.
cfg_def.lfp_chans = []; % lfp channels to be plotted can be empty. Can be 1 or more, but best to keep it less than 3. should be rows in csc.data.
cfg_def.resize = 1; % use the gui to select the data for resizing.
cfg_def.state_keys = {'rightarrow','uparrow', 'downarrow', 'leftarrow', 'numpad0', 'backspace','backquote' }; % which key to press for each state_val
cfg_def.state_name = {'Wake', 'SWS', 'REM', 'Quiescence', 'Transition', 'Redo', 'Exit'}; %
cfg_def.state_val = 1:length(cfg_def.state_name); % what numerical value for each state.
cfg_def.spec.win_s = 2^9; % spectrogram window size.
cfg_def.spec.onverlap = pow2(floor(log2(cfg_def.spec.win_s/4)));%cfg_def.spec.win_s /8 ; % overlap
cfg_def.spec.freq_low = .5:0.25:14; % frequency range for spectrogram for delta/theta
cfg_def.spec.freq_high = 30:5:250; % frequency range for spectrogram in the SWR range
cfg_def.spec.lfp_chan = 1; % which channel to use for the spectrogram
cfg_def.smooth_spec = 2; 
cfg_def.Fs = mode(diff(tvec)); % best guess of the sampling frequency.
cfg_def.delta_f = [1 4];
cfg_def.theta_f = [6 9];


cfg = ProcessConfig(cfg_def, cfg_in);

%% setup some variables
if isempty(score_in)
    score = NaN(size(tvec));
else
    score = score_in;
end


%% generate spectrograms

for iChan = length(cfg.lfp_chans):-1:1
    fprintf('<strong>%s</strong>: processing data channel %d\n', mfilename, iChan)
    [~, F{iChan}, T{iChan},P{iChan}] =  spectrogram(data(cfg.lfp_chans(iChan),:),hanning(cfg_def.spec.win_s),cfg.spec.onverlap,cfg.spec.freq_low,1/cfg.Fs);
    

    
%     if cfg.smooth_spec
%         P{iChan} = smoothdata(P{iChan}, 2, 'gaussian',  2);
%     end
    
    fprintf('<strong>%s</strong>: processing data channel %d\n', mfilename, iChan)
    [~, Fh{iChan}, Th{iChan},Ph{iChan}] =  spectrogram(data(cfg.lfp_chans(iChan),:),hanning(cfg_def.spec.win_s/4),cfg.spec.onverlap/8,cfg.spec.freq_high,1/cfg.Fs);
%     if cfg.smooth_spec
%         Ph{iChan} = smoothdata(Ph{iChan}, 2, 'gaussian',  2);
%     end

% remove artifacts or SVSs
%     art = zscore(abs(sum(P{iChan})));
%     art_idx = art>4; 
%     art = zscore(abs(sum(Ph{iChan})));
%     art_idx_h  = art>4; 
%     all_art = art_idx==1 & art_idx_h==1; 
%     
%     Ph{iChan}(:,all_art) = NaN;
%     P{iChan}(:,all_art) = NaN; 
end

% % buz method
%    freqlist = logspace(0,2,100);
%     [swFFTspec,swFFTfreqs,t_clus] = spectrogram(data(cfg.lfp_chans(iChan),:),512,256,freqlist,1/cfg.Fs);
%      swFFTspec = abs(swFFTspec);
%     specdt = mode(diff(t_clus));
%     zFFTspec = zscore(log10(swFFTspec)');
%     % Remove transients before calculating SW histogram
%     totz = zscore(abs(sum(zFFTspec')));
%     badtimes = find(totz>5);
% %     swFFTspec(badtimes,:) = NaN; 
% imagesc(t_clus,swFFTfreqs,swFFTspec)
%  axis xy
%         set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
%         set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
% % axis xy
%% filtered signals? [maybe add?]
% % delta filter.
% cfg_filt_d = [];
% cfg_filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
% cfg_filt_d.f  = cfg.delta_f;
% cfg_filt_d.order = 8; %type filter order
% % cfg_filt_d.display_filter = 1; % use this to see the fvtool and wait for a keypress before continuing.
% delta_csc = FilterLFP(cfg_filt_d,tsd(tvec, data));
% close all
%
%
% % filter into the theta band
% cfg_filt_t = [];
% cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
% cfg_filt_t.f  = cfg.theta_f;
% cfg_filt_t.order = 3; %type filter order
% % cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% % reason.  .
% theta_csc = FilterLFP(cfg_filt_t, lfp);

%% create basic plot with channel(s) of interest + emg if available.
hFig = figure;

if ~isempty(score_in)
    Spec_subs = 2:length(cfg.lfp_chans):((length(cfg.lfp_chans)*2)*(length(cfg.lfp_chans)*2)+1);
else
    Spec_subs = 1:length(cfg.lfp_chans):((length(cfg.lfp_chans)*2)*(length(cfg.lfp_chans)*2));
end
% add another subplot if you have an emg.

if ~isempty(emg)
    nSubplots = Spec_subs(end)+2;
else
    nSubplots = Spec_subs(end)+1;
end

if ~isempty(score_in)
    ax(1000) = subplot(nSubplots, 1, 1);
    plot(tvec- tvec(1), score_in);
    text(tvec(1)- tvec(1), median(score_in), 'Score')
end
% plot an empty sleep score. [to do] make this fill in with the score_in
if ~isempty(score_in)
    ax(1000) = subplot(nSubplots, 1, 1);
    plot(tvec- tvec(1), score_in);
    text(tvec(1)- tvec(1), median(score_in), 'Score')
end


for iSub = 1:length(cfg.lfp_chans)
    ax(iSub) = subplot(nSubplots, 1, Spec_subs(iSub):Spec_subs(iSub+1));
    h = imagesc(Th{iSub}, Fh{iSub}, 10*log10(Ph{iSub}));
    set(gca, 'xtick', [])
    axis xy
    
    ax(iSub+100) = subplot(nSubplots, 1, Spec_subs(iSub)+2:Spec_subs(iSub)+3);
    h = imagesc(T{iSub}, F{iSub}, 10*log10(P{iSub}));
    ylim([F{iSub}(1) F{iSub}(end)])
    axis xy
    
end
ax(99) = subplot(nSubplots, 1, nSubplots-1);
plot(tvec- tvec(1), data(cfg.lfp_chans(iChan),:))
%     xlim(cfg.emg_range)
text(tvec(1)- tvec(1), median(data(cfg.lfp_chans(iChan),:)), 'LFP')

if ~isempty(emg)
    ax(999) = subplot(nSubplots, 1, nSubplots);
    plot(tvec- tvec(1), emg(cfg.emg_chan,:))
%     ylim([min(emg(cfg.emg_chan,:)), max(emg(cfg.emg_chan,:))])
    ylim(cfg.emg_range)
    text(tvec(1)- tvec(1), cfg.emg_range(2)*.8, 'EMG')
end

linkaxes(ax, 'x')
xlim([T{1}(1) T{1}(end)])

[cut_x,~] = ginput(1);
if cut_x < 0
    cut_x = 0;
end
%% GUI

% ButtonH=uicontrol('Parent',hFig,'Style','pushbutton','String','Next','Units','normalized','Position',[0.0 0.5 0.4 0.2],'Visible','on');

%% start the scrolling and scoring
% input('Paused for inspection. press any key to continue')

% actively plot the data. 
hFig2 = figure;
    this_idx = nearest_idx([cut_x cut_x+cfg.tvec_range(2)],tvec-tvec(1)); % get the tvec index for the current window
    spec_idx = nearest_idx([cut_x cut_x+cfg.tvec_range(2)], T{1});
    spec_h_idx = nearest_idx([cut_x cut_x+cfg.tvec_range(2)], Th{1});

% plot an empty sleep score. [to do] make this fill in with the score_in
if ~isempty(score_in)
    ax(1000) = subplot(nSubplots, 1, 1);
    plot(tvec(this_idx(1):this_idx(2))- tvec(1), score_in(this_idx(1):this_idx(2)));
%     text(tvec(1)- tvec(1), median(score_in), 'Score')
    xlim([min(score_in), max(score_in)]);
    drawnow;
end

for iSub = 1:length(cfg.lfp_chans)
    ax(iSub) = subplot(nSubplots, 1, Spec_subs(iSub):Spec_subs(iSub+1));
    h = imagesc(Th{iSub}(spec_h_idx(1):spec_h_idx(2)), Fh{iSub}, 10*log10(Ph{iSub}(:,spec_h_idx(1):spec_h_idx(2))));
    set(gca, 'xtick', [])
    axis xy
    drawnow;
    
    ax(iSub+100) = subplot(nSubplots, 1, Spec_subs(iSub)+2:Spec_subs(iSub)+3);
    h = imagesc(T{iSub}(spec_idx(1):spec_idx(2)), F{iSub}, 10*log10(P{iSub}(:,spec_idx(1):spec_idx(2))));
%     ylim([F{iSub}(1) ceil(F{iSub}(end)*1.1)])
    axis xy
    drawnow;
end
ax(99) = subplot(nSubplots, 1, nSubplots-1);
plot(tvec(this_idx(1):this_idx(2))- tvec(1), data(cfg.lfp_chans(iChan),this_idx(1):this_idx(2)))
%     xlim(cfg.emg_range)
% text(tvec(1)- tvec(1), median(data(cfg.lfp_chans(iChan),:)), 'LFP')
drawnow; 

if ~isempty(emg)
    ax(999) = subplot(nSubplots, 1, nSubplots);
    plot(tvec(this_idx(1):this_idx(2)) - tvec(1), emg(cfg.emg_chan,this_idx(1):this_idx(2)))
        ylim(cfg.emg_range)
%     ylim([min(emg(cfg.emg_chan,:)), max(emg(cfg.emg_chan,:))])
%     text(tvec(1)- tvec(1), cfg.emg_range(2)*.8, 'EMG')
drawnow; 
end

linkaxes(ax, 'x')
% xlim([T{1}(1) T{1}(end)])

for iL = 1:length(cfg.state_keys)
    fprintf('%s  = %s    \n', cfg.state_name{iL}, cfg.state_keys{iL})
end


xlim([cut_x cut_x+cfg.tvec_range(2)])
xval = xlim;


done =0; 
while ishandle(h)
    
    % to do make this more specific with selectable buttons and
    % roi = drawline('StripeColor','r');
    % if exist(
    % start = roi.Position(1,1);
    % stop = roi.Position(2,1);
    % this_idx = nearest_idx([start stop],tvec-tvec(1)); % get the tvec index for the current window
%     this_idx = nearest_idx(xlim,tvec-tvec(1)); % get the tvec index for the current window
    
    
    was_a_key = waitforbuttonpress;
    key_hit = get(gcf, 'CurrentKey');
    if ismember(key_hit, cfg.state_keys(1:end-2))
        iK = find(ismember(cfg.state_keys, key_hit));
        fprintf('%s ', cfg.state_name{iK})
        score(this_idx(1):this_idx(2)) = cfg.state_val(iK);
        x_lim = xlim;
        fprintf('     %.0f   %.0f    %0.2f%%\n', x_lim(1), x_lim(2), round(x_lim(2)/(tvec(end)-tvec(1)),2)*100)
        if done == 1 
            disp('Exiting')
            break
        end
        %         subplot(nSubplots, 1, 1)
        %         plot(tvec- tvec(1), score);
        %         drawnow;
        if x_lim(2)+cfg.tvec_range(2) < T{1}(end)
            this_idx = nearest_idx(xlim+cfg.tvec_range(2),tvec-tvec(1)); % get the tvec index for the current window
            spec_idx = nearest_idx(xlim+cfg.tvec_range(2), T{1});
            spec_h_idx = nearest_idx(xlim+cfg.tvec_range(2), Th{1});
        else
            disp('at the end')
            done = 1; % used to break loop. 
            this_idx = [nearest_idx(x_lim(1)+cfg.tvec_range(2),tvec-tvec(1)), length(tvec)]; % get the tvec index for the current window
            spec_idx = [nearest_idx(x_lim(1)+cfg.tvec_range(2), T{1}),length(T{1})];
            spec_h_idx = [nearest_idx(x_lim(1)+cfg.tvec_range(2), Th{1}),length(Th{1})];
        end
    elseif ismember(key_hit,  cfg.state_keys(end))
        disp('Exiting')
        break
    elseif ismember(key_hit,  cfg.state_keys(end-1))
                fprintf('REDO   \n')
        this_idx = nearest_idx(xlim-cfg.tvec_range(2),tvec-tvec(1)); % get the tvec index for the current window
        spec_idx = nearest_idx(xlim-cfg.tvec_range(2), T{1});
        spec_h_idx = nearest_idx(xlim-cfg.tvec_range(2), Th{1});
    end
    


% plot an empty sleep score. [to do] make this fill in with the score_in
if ~isempty(score_in)
    ax(1000) = subplot(nSubplots, 1, 1);
    plot(tvec(this_idx(1):this_idx(2))- tvec(1), score_in(this_idx(1):this_idx(2)));
%     text(tvec(1)- tvec(1), median(score_in), 'Score')
    xlim([min(score_in), max(score_in)]);
    drawnow;
end

for iSub = 1:length(cfg.lfp_chans)
    ax(iSub) = subplot(nSubplots, 1, Spec_subs(iSub):Spec_subs(iSub+1));
    h = imagesc(Th{iSub}(spec_h_idx(1):spec_h_idx(2)), Fh{iSub}, 10*log10(Ph{iSub}(:,spec_h_idx(1):spec_h_idx(2))));
    set(gca, 'xtick', [])
    axis xy
    drawnow;
    
    ax(iSub+100) = subplot(nSubplots, 1, Spec_subs(iSub)+2:Spec_subs(iSub)+3);
    h = imagesc(T{iSub}(spec_idx(1):spec_idx(2)), F{iSub}, 10*log10(P{iSub}(:,spec_idx(1):spec_idx(2))));
%     ylim([F{iSub}(1) ceil(F{iSub}(end)*1.1)])
    axis xy
    drawnow;
end
ax(99) = subplot(nSubplots, 1, nSubplots-1);
plot(tvec(this_idx(1):this_idx(2))- tvec(1), data(cfg.lfp_chans(iChan),this_idx(1):this_idx(2)))
xlim([tvec(this_idx(1))- tvec(1), tvec(this_idx(2))- tvec(1)]);
% text(tvec(1)- tvec(1), median(data(cfg.lfp_chans(iChan),:)), 'LFP')
drawnow; 

if ~isempty(emg)
    ax(999) = subplot(nSubplots, 1, nSubplots);
    plot(tvec(this_idx(1):this_idx(2)) - tvec(1), emg(cfg.emg_chan,this_idx(1):this_idx(2)))
    ylim(cfg.emg_range)
    xlim([tvec(this_idx(1))- tvec(1), tvec(this_idx(2))- tvec(1)]);
%     text(tvec(1)- tvec(1), cfg.emg_range(2)*.8, 'EMG')
drawnow; 
end
%     xval = xval+cfg.tvec_range(2);
%     xlim(xval)
%     if ~isempty(emg)
%         ax(999) = subplot(nSubplots, 1, nSubplots);
%         ylim([min(emg(cfg.emg_chan,:)), max(emg(cfg.emg_chan,:))])
%     end
end

%% [to do] get the transitions and convert to iv format. 

fprintf('<strong>%s</strong>: complete. Returning complete sleep score\n', mfilename); 




