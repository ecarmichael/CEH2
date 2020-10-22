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


% conver the data to the ft format
data_ft = MS_TSDtoFT([], csc); % convert to ft format.

% convert to trials
cfg_trl = [];
cfg_trl.t = cat(1,events);
cfg_trl.t = cfg_trl.t - data_ft.hdr.FirstTimeStamp;
cfg_trl.twin = [-3 5];
cfg_trl.hdr = data_ft.hdr;

trl = ft_maketrl(cfg_trl);

cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data_ft);


%% setup the initial analysis.

cfg              = []; % start with empty cfg
cfg.output       = 'pow';
cfg.channel      = data_ft.label{1};
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 20:.2:80; % frequencies of interest
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;%20./cfg.foi;  % window size: fixed at 0.5s
cfg.toi          = -2.5:0.05:4.5; % times of interest
%cfg.pad          = 'nextpow2'; % recommened by FT to make FFT more efficient.

TFR = ft_freqanalysis(cfg, data_trl);

%% track config for plotting.
fprintf('Spec using %0.0d swrs. Method: %s, Taper: %s', length(trl),cfg.method, cfg.taper);
if ~isempty(label)
freq_params_str = sprintf(label);
else
    freq_params_str = sprintf('Spec using %0.0d swrs. Method: %s, Taper: %s', length(trl),cfg.method, cfg.taper);
end
cfg = [];
cfg.channel      = data_ft.label{1};
cfg.baseline     = [-2 -.01];
cfg.baselinetype = 'relative';
cfg.title = freq_params_str;
ft_singleplotTFR(cfg, TFR);

vline(0)