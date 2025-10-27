function [hypno, csc, emg] = MS_get_hypno(csc, emg, wake_idx, emg_rms_prctile, d_t_ratio)
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

if nargin < 1
    Meta = MS_Load_meta;
    cfg_csc.fc = {Meta.goodCSC};%, 'CSC3.ncs'};  % pick the channels to load
    csc = LoadCSC(cfg_csc); % load the csc data
    
    cfg_csc.fc = {Meta.EMG};%, 'CSC3.ncs'};  % pick the channels to load
    emg = LoadCSC(cfg_csc); % load the csc data
    emg_rms_prctile = 90;
    d_t_ratio = 1;

elseif nargin < 2
    Meta = MS_Load_meta;
    cfg_csc.fc = {Meta.EMG};%, 'CSC3.ncs'};  % pick the channels to load
    emg = LoadCSC(cfg_csc); % load the csc data
    emg_rms_prctile = 90;
    d_t_ratio = 1;

elseif nargin < 3
    wake_idx = [];
    emg_rms_prctile = 90; 
        d_t_ratio = 1;

elseif nargin < 4
    emg_rms_prctile = 90; 
        d_t_ratio = 1;
elseif nargin < 5
        d_t_ratio = 1;
end

if sum(diff(unique(diff(csc.tvec)))>csc.cfg.hdr{1}.SamplingFrequency*4) > 1
    error('Data contains gaps. Should be continuous')
end

cord = linspecer(5);

% d_t_ratio = 1.5;
% emg_rms_prctile = 90;




%% mark out wake periods from the EMG;
%This
sat_idx = false(1, length(csc.data));%(abs(csc.data) >= prctile(abs(csc.data), 90));% | (csc.data <= prctile(csc.data, 5));

if ~isempty(wake_idx)
    for ii = 1:size(wake_idx,1)
        sat_idx(wake_idx(ii, 1):wake_idx(ii, 2)) = 1;
    end
end




%% get the EMG RMS
if isstruct(emg)
    emg_rms = sqrt(movmedian(emg.data.^2, csc.cfg.hdr{1}.SamplingFrequency*10));      % RMS Value Over ‘WinLen’ Samples
    emg_og = emg_rms;
else
%     emg_rms_prctile = 90;
    emg_rms = movmedian(emg, csc.cfg.hdr{1}.SamplingFrequency*10) ;
    emg_og = emg_rms;
end

emg_rms(sat_idx) = NaN;
% emg_rms_z = zscore(emg_rms);

% emg_rms = movmedian(emg_rms,(csc.cfg.hdr{1}.SamplingFrequency)*5);

move_thresh = prctile(emg_rms, emg_rms_prctile);
% REM_thresh = prctile(emg_rms, 15);

move_idx = emg_rms >move_thresh;
%% get the theta delta ratio
cfg_con_t = [];
cfg_con_t.threshold = 0;
cfg_con_t.f = [6 12]; cfg_con_t.type = 'fdesign';
%     cfg_con_t.display_filter = 1
csc_th = FilterLFP(cfg_con_t,csc);
% csc_th.data = smooth(csc_th.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_th.data = smoothdata(abs(hilbert(csc_th.data)), 'movmean', csc.cfg.hdr{1}.SamplingFrequency*10);

cfg_con_d = [];
cfg_con_d.threshold = 0;
cfg_con_d.f = [1.5 4]; cfg_con_d.type = 'fdesign';
% cfg_con_d.display_filter = 1
csc_delta = FilterLFP(cfg_con_d,csc);
% csc_delta.data = smooth(csc_delta.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_delta.data = smoothdata(abs(hilbert(csc_delta.data)), 'movmean', csc.cfg.hdr{1}.SamplingFrequency*10);


% block out any known wake/problematic periods
% csc_th.data(sat_idx) = NaN;
% csc_delta.data(sat_idx) = NaN;

ratio = csc;
ratio.data = zscore((csc_th.data ./ csc_delta.data)');

z_ratio = nan(1, length(csc.tvec));
z_ratio(~move_idx|sat_idx) = zscore((csc_th.data(1, ~move_idx|sat_idx)./emg_rms(1,~move_idx|sat_idx)));%csc_delta.data(~move_idx)));

z_ratio = ratio.data';

hypno_init = ones(size(csc.tvec))*2; % set everything to SWS;
hypno_init(~move_idx & (z_ratio > d_t_ratio)) = 3; % putative REM
hypno_init(move_idx | sat_idx) = 1; %awake based on movement;

labels = {'Wake', 'SWS', 'REM'};

%% smooth things out to avoid jumps.

% sws
wake_tsd = tsd(csc.tvec, (hypno_init == 1)', 'wake');

wake_cfg.threshold = 0;
wake_cfg.dcn = '>';
wake_cfg.minlen = 5;
wake_cfg.merge_thr = 5;
wake_iv = TSDtoIV(wake_cfg, wake_tsd);


REM_tsd = tsd(csc.tvec, (hypno_init == 3)', 'REM');

rem_cfg.threshold = 0;
rem_cfg.dcn = '>';
rem_cfg.minlen = 10;
rem_cfg.merge_thr = 5;
REM_iv = TSDtoIV(rem_cfg, REM_tsd);

SWS_tsd = tsd(csc.tvec, (hypno_init == 2)', 'SWS');

sws_cfg.threshold = 0;
sws_cfg.dcn = '>';
sws_cfg.minlen = 30;
sws_cfg.merge_thr = 5;
SWS_iv = TSDtoIV(sws_cfg, SWS_tsd);

% cfg = [];
% cfg.display = 'iv';
% PlotTSDfromIV(cfg, REM_iv, csc)

% rebuild hypno_init with cut REM

hypno_out = ones(size(csc.tvec))*2; % set everything to SWS;

for ii = 1:length(wake_iv.tstart)
    hypno_out(nearest_idx3(wake_iv.tstart(ii), csc.tvec):  nearest_idx3(wake_iv.tend(ii), csc.tvec)) = 1;
end


for ii = 1:length(REM_iv.tstart)
    hypno_out(nearest_idx3(REM_iv.tstart(ii), csc.tvec):  nearest_idx3(REM_iv.tend(ii), csc.tvec)) = 3;
end

% hypno_init(move_idx) = 1; %awake based on movement;



%% check figure

% decimate data for speed
deci_fact = 10; 


        csc_plot.data = decimate(csc.data,deci_fact);
        csc_plot.tvec = csc.tvec(1:deci_fact:end);
        csc_plot.cfg.hdr{1}.SamplingFrequency = csc.cfg.hdr{1}.SamplingFrequency./deci_fact;

        emg_plot = emg_rms(1:deci_fact:end); % subsamples instead of decimate due to NaNs from saturations. 
        move_plot = move_idx(1:deci_fact:end); 
        z_ratio_plot= decimate(z_ratio,deci_fact);
        ratio_plot.data = decimate(ratio.data,deci_fact);
        hypno_out_plot = hypno_out(1:deci_fact:end); 
        
figure
clf
ax(1) = subplot(7,1,1:4);
cla
yyaxis right
hold on
% plot((csc.tvec - csc.tvec(1)),  emg_og'*1000,  'color', [0.3 0.3 0.3]);
plot((csc_plot.tvec - csc_plot.tvec(1)),  emg_plot'*1000,  'color', cord(2,:));

ylim([min(emg_plot*1000), 3*max(emg_plot*1000)])
yline(move_thresh*1000)
set(gca, 'xtick', [])
ylabel('EMG voltage (mV)')

% ylim([0 50])
yyaxis left
plot((csc_plot.tvec - csc_plot.tvec(1)), csc_plot.data*1000, 'color', cord(1,:));
legend({'Raw LFP', 'smooth EMG-rms'});
% xlim([min((csc.tvec - csc.tvec(1))/60/60) max((csc.tvec - csc.tvec(1))/60/60)])
ylim([.75*min(csc_plot.data*1000), 1.25*max(csc_plot.data*1000)])
ylabel('LFP voltage')

set(gca, 'xtick', [])
% ylabel('LFP voltage')
title(['Auto-Hypno (note: decimated by ' num2str(deci_fact) ')'])



ax(2) = subplot(7,1,5:6);
cla
hold on
plot((csc_plot.tvec - csc_plot.tvec(1)),  ~move_plot & (z_ratio_plot > .5), 'k');
plot((csc_plot.tvec - csc_plot.tvec(1)), ratio_plot.data,  'color', cord(3,:));
% plot((csc.tvec - csc.tvec(1)), sat_idx - 2,  'color', cord(5,:));
ylabel('Theta-delta ratio')
legend({'REM idx', 'theta/delta z'});

% yline(.5);
set(gca, 'xtick', [])

ax(3) = subplot(7,1,7);
cla
hold on
vec = zeros(size(hypno_out_plot));
plot((csc_plot.tvec(hypno_out_plot == 2) - csc_plot.tvec(1)), vec(hypno_out_plot == 2)+2,'s','MarkerEdgeColor', cord(2,:),'MarkerFaceColor',cord(2,:), 'linewidth', 1)
plot((csc_plot.tvec(hypno_out_plot == 3) - csc_plot.tvec(1)), vec(hypno_out_plot == 3)+3,'s','MarkerEdgeColor', cord(3,:),'MarkerFaceColor',cord(3,:), 'linewidth', 1)
plot((csc_plot.tvec(hypno_out_plot == 1) - csc_plot.tvec(1)), vec(hypno_out_plot == 1)+1,'s','MarkerEdgeColor', cord(1,:),'MarkerFaceColor',cord(1,:), 'linewidth', 1)

% imagesc((csc.tvec - csc.tvec(1)), 1, hypno_init')
ytickValues = [1, 2, 3];  % Specify the values where you want tick marks
yticks(ytickValues);
set(gca, 'YTickLabel', labels)
% colormap(linspecer(3));
% cb=colorbar;
% cb.Position(1) = cb.Position(1) + 1e-1;
% cb.Ticks = [1:3];
% cb.TickLabels = {'Wake', 'SWS', 'REM'};
xlabel('time (s)')

% legend('SWS', 'REM', 'Awake')

linkaxes(ax, 'x')
xlim([min((csc_plot.tvec - csc_plot.tvec(1))) max((csc_plot.tvec - csc_plot.tvec(1)))])

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
cfg_hypno.SWS = SWS_iv;


hypno = [];
hypno.data = hypno_out;
hypno.tvec = csc.tvec;
hypno.labels = labels;
hypno.cfg = cfg_hypno;
