function Triggered_Spec_FT(csc, events, label)
%% Triggered_Spec_FT: wrapper function for generating an event triggered spectrogram using the FieldTrip ft_freqanalysis.
%
%
%
%    Inputs:
%    - csc [struct]  in the
%
%    - events [1 x nEvent start times
%
%    Outputs:
%    - none atm.
%
%
%
%
% EC 2020-10-20   initial version.  no flexibility.
%
% TODO: add cfg configs
%
%% initialize
if nargin == 2
    label = [];
elseif nargin <2
    fprintf('<strong>%s<\strong> Requires csc and events inputs', mfilename);
end

if max(diff(csc.tvec)) > csc.cfg.hdr{1}.SamplingFrequency
    error('discontinous data this may cause problems')
end
% conver the data to the ft format
data_ft = MS_TSDtoFT([], csc); % convert to ft format.

% convert to trials
cfg_trl = [];
cfg_trl.t = cat(1,events);
cfg_trl.t = cfg_trl.t - data_ft.hdr.FirstTimeStamp;
cfg_trl.twin = [-5 5];
cfg_trl.hdr = data_ft.hdr;

trl = ft_maketrl(cfg_trl);

cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data_ft);


%% setup the initial analysis.

cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = data_ft.label{1};
cfg.method       = 'wavelet';
cfg.taper        = 'hanning';
cfg.foi          = 2:.5:120; % frequencies of interest
% cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;%20./cfg.foi;  % window size: fixed at 0.5s
cfg.toi          = -5:((1/data_trl.fsample)*10):3; % times of interest
cfg.pad          = 'nextpow2'; % recommened by FT to make FFT more efficient.

TFR = ft_freqanalysis(cfg, data_trl);

%% track config for plotting.
fprintf('Spec using %0.0d events. Method: %s', length(trl),cfg.method);
if ~isempty(label)
freq_params_str = sprintf(label);
else
    freq_params_str = sprintf('Spec using %0.0d swrs. Method: %s', length(trl),cfg.method);
end
cfg = [];
cfg.channel      = data_ft.label{1};
cfg.baseline     = [-5 -2];
cfg.baselinetype = 'zscore';
cfg.title = freq_params_str;
ft_singleplotTFR(cfg, TFR);

vline(0, '--k')