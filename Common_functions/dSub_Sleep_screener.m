function [hypno, csc, emg] = dSub_Sleep_screener(csc, emg)
%% dSub_Sleep_screener: run an intial screening of sleep ephys data to classify data as movement or rough sleep stages. 
%
%
%
%    Inputs: 
%    - csc: [struct] LFP data for the channel of interest. in the TSD

%    format from MS_LoadCSC; 
%    
%     - emg: [struct] emg data for the channel of interest. in the TSD
%    format from MS_LoadCSC; 
%
%    Outputs: 
%    - hypno_init: [struct]  
%           data = [1 x nSamples] initial sleep scoring. 
%           labels: [cells] labels for the data values. should be 1= wake,
%               2 = SWS, 3 = REM
%           cfg: [stgruct] contains configs and IVs; 
%
%    - csc: [struct] LFP data used. 
%
%    - emg: [struct] EMG used. 
%
%
% EC 2022-08-22   initial version 
%
%% initialize
    Meta = MS_Load_meta;

if nargin < 1
    cfg_csc.fc = {Meta.goodCSC};%, 'CSC3.ncs'};  % pick the channels to load
    csc = LoadCSC(cfg_csc); % load the csc data
    
    cfg_csc.fc = {Meta.EMG};%, 'CSC3.ncs'};  % pick the channels to load
    emg = LoadCSC(cfg_csc); % load the csc data
elseif nargin < 2
    
    cfg_csc.fc = {Meta.EMG};%, 'CSC3.ncs'};  % pick the channels to load
    emg = LoadCSC(cfg_csc); % load the csc data
end

if sum(diff(unique(diff(csc.tvec)))>csc.cfg.hdr{1}.SamplingFrequency*4) > 1
    error('Data contains gaps. Should be continuous')
end

cord = linspecer(5);

d_t_ratio = .5; 
emg_rms_prctile = 70; 
%% get the EMG RMS

emg_rms = sqrt(movmedian(emg.data.^2, csc.cfg.hdr{1}.SamplingFrequency*5));      % RMS Value Over ‘WinLen’ Samples

% emg_rms_z = zscore(emg_rms); 

% emg_rms = movmedian(emg_rms,(csc.cfg.hdr{1}.SamplingFrequency)*5);

move_thresh = prctile(emg_rms, emg_rms_prctile);
% REM_thresh = prctile(emg_rms, 15); 

move_idx = emg_rms >move_thresh; 
%% get the theta delta ratio
cfg_con_t = [];
cfg_con_t.threshold = 0;
cfg_con_t.f = [6 12]; cfg_con_t.type = 'fdesign';
%     cfg_con_f.display_filter = 1
csc_th = FilterLFP(cfg_con_t,csc);
% csc_th.data = smooth(csc_th.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_th.data = smooth(abs(hilbert(csc_th.data)), csc.cfg.hdr{1}.SamplingFrequency*10);

cfg_con_d = [];
cfg_con_d.threshold = 0;
cfg_con_d.f = [1 4]; cfg_con_d.type = 'fdesign';
%     cfg_con_f.display_filter = 1
csc_delta = FilterLFP(cfg_con_d,csc);
% csc_delta.data = smooth(csc_delta.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_delta.data = smooth(abs(hilbert(csc_delta.data)), csc.cfg.hdr{1}.SamplingFrequency*10);


ratio = csc;
ratio.data = zscore((csc_th.data ./ csc_delta.data)');

z_ratio = nan(1, length(csc.tvec));
z_ratio(~move_idx) = zscore((csc_th.data(~move_idx)./emg_rms(~move_idx)'));%csc_delta.data(~move_idx)));


sat_idx = (csc.data == max(csc.data)) | (csc.data == min(csc.data));

hypno_init = ones(size(csc.tvec))*2; % set everything to SWS; 
hypno_init(~move_idx & z_ratio > d_t_ratio) = 3; % putative REM
hypno_init(move_idx) = 1; %awake based on movement; 

labels = {'Wake', 'SWS', 'REM'}; 

%% smooth things out to avoid jumps. 

% sws
wake_tsd = tsd(csc.tvec, (hypno_init == 1)', 'wake');

wake_cfg.threshold = 0; 
wake_cfg.dcn = '>';
wake_cfg.minlen = 10;
wake_cfg.merge_thr = 10; 
wake_iv = TSDtoIV(wake_cfg, wake_tsd);


REM_tsd = tsd(csc.tvec, (hypno_init == 3)', 'REM');

rem_cfg.threshold = 0; 
rem_cfg.dcn = '>';
rem_cfg.minlen = 15;
rem_cfg.merge_thr = 1; 
REM_iv = TSDtoIV(rem_cfg, REM_tsd);
 
% cfg = [];
% cfg.display = 'iv';
% PlotTSDfromIV(cfg, REM_iv, csc)

% rebuild hypno_init with cut REM

hypno_out = ones(size(csc.tvec))*2; % set everything to wake; 

for ii = 1:length(wake_iv.tstart)
   hypno_out(nearest_idx3(wake_iv.tstart(ii), csc.tvec):  nearest_idx3(wake_iv.tend(ii), csc.tvec)) = 1; 
end


for ii = 1:length(REM_iv.tstart)
   hypno_out(nearest_idx3(REM_iv.tstart(ii), csc.tvec):  nearest_idx3(REM_iv.tend(ii), csc.tvec)) = 3; 
end

% hypno_init(move_idx) = 1; %awake based on movement; 



%% check figure
figure(222)

ax(1) = subplot(7,1,1:4);
cla
tic
hold on
yyaxis right
plot((csc.tvec - csc.tvec(1)),  emg_rms*1000,  'color', cord(2,:));
ylim([min(emg_rms*1000), 3*max(emg_rms*1000)])
yline(move_thresh*1000)
set(gca, 'xtick', [])

% ylim([0 50])
yyaxis left
plot((csc.tvec - csc.tvec(1)), csc.data*1000, 'color', cord(1,:));
legend({'Raw LFP', 'smooth EMG-rms'});
% xlim([min((csc.tvec - csc.tvec(1))/60/60) max((csc.tvec - csc.tvec(1))/60/60)])
ylim([.75*min(csc.data*1000), 1.25*max(csc.data*1000)])

set(gca, 'xtick', [])
% ylabel('LFP voltage')
toc

ax(2) = subplot(7,1,5:6);
cla
hold on
plot((csc.tvec - csc.tvec(1)),  ~move_idx & z_ratio > .5, 'k');
plot((csc.tvec - csc.tvec(1)), z_ratio,  'color', cord(3,:));
legend({'excluded idx', 'theta/delta z'});

yline(.5); 
set(gca, 'xtick', [])

ax(3) = subplot(7,1,7);
cla
hold on
vec = zeros(size(hypno_out)); 
plot((csc.tvec(hypno_out == 2) - csc.tvec(1)), vec(hypno_out == 2)+2,'s','MarkerEdgeColor', cord(2,:),'MarkerFaceColor',cord(2,:), 'linewidth', 1)
plot((csc.tvec(hypno_out == 3) - csc.tvec(1)), vec(hypno_out == 3)+3,'s','MarkerEdgeColor', cord(3,:),'MarkerFaceColor',cord(3,:), 'linewidth', 1)
plot((csc.tvec(hypno_out == 1) - csc.tvec(1)), vec(hypno_out == 1)+1,'s','MarkerEdgeColor', cord(1,:),'MarkerFaceColor',cord(1,:), 'linewidth', 1)

    % imagesc((csc.tvec - csc.tvec(1)), 1, hypno_init')
    set(gca, 'YTickLabel', {'Wake', 'SWS', 'REM'})
% colormap(linspecer(3)); 
% cb=colorbar;
% cb.Position(1) = cb.Position(1) + 1e-1;
% cb.Ticks = [1:3]; 
% cb.TickLabels = {'Wake', 'SWS', 'REM'}; 
xlabel('time (s)')

% legend('SWS', 'REM', 'Awake')

linkaxes(ax, 'x')
xlim([min((csc.tvec - csc.tvec(1))) max((csc.tvec - csc.tvec(1)))])

% SetFigure([], gcf);
maximize
%% collect the values

cfg_hypno = [];
cfg_hypno.d_t_ratio = d_t_ratio; 
cfg_hypno.emg_rms_prctile = emg_rms_prctile; 
cfg_hypno.cfg_con_d = cfg_con_d; 
cfg_hypno.cfg_con_t = cfg_con_t; 
cfg_hypno.wake = wake_iv; 
cfg_hypno.REM = REM_iv; 

hypno = [];
hypno.data = hypno_out; 
hypno.tvec = csc.tvec; 
hypno.labels = labels;
hypno.cfg = cfg_hypno; 