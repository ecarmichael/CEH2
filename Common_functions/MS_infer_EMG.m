function corr_out = MS_infer_EMG(cfg_in, LFP_1, LFP_2)
%% MS_infer_EMG: infer the EMG signal from high frequency xcorr between two LFP channels (in the same structure is probably best). Based on the methods from the Buzsaki lab ('bz_EMGFromLFP')
%
%
%
%    Inputs: 
%    - cfg_in  [struct] configuration parameters to overwrite the defaults
%    below.
%
%    - LFP1  [struct] TSD LFP structure from LoadCSC to be correlated with
%    LFP2
%
%    - LFP2 optional [struct] TSD LFP structure from LoadCSC to be correlated with
%    LFP1
%
%
%
%    Outputs: 
%    - EMG_est:  [struct]   EMG estimate in the TSD format. 
%
%
% This method is based on 'bz_EMGFromLFP' by Erik Schomburg and others
% (Watson, Levenstein, Tingley & Swanson) from 
% https://github.com/buzsakilab/buzcode/blob/05a8bd55ac48a6d50069c3f392d55fdc6e8cd5ec/detectors/bz_EMGFromLFP.m
%
% EC 2020-11-13   initial version to fit the LFP TSD format. 
%
%
%
%% initialize

if nargin == 2 % if input is one TSD then pull out the two signals.
    if size(LFP_1.data)  ==1
        error('If only using 2 inputs, then LFP_1 needs at least 2 data channels.')
    end
    % rebuild the LFP_2 into its own TSD structure (not efficient but
    % needed for filtering. 
    LFP_2 = LFP_1; 
    LFP_2.data = [];
    LFP_2.data = LFP_1.data(2,:);
    LFP_2.label = LFP_1.label{2}; 
    LFP_2.cfg.hdr = LFP_1.cfg.hdr(2);
    % remove LFP_2 from the original LFP_1
    LFP_1.data = LFP_1.data(1,:); % take out the 
        LFP_1.label = LFP_1.label{1}; 
    LFP_1.cfg.hdr = LFP_1.cfg.hdr(1);
end


cfg_def = [];
cfg_def.freq_bands = [300 600]; % bandpass 
cfg_def.window = [-.5 .5]; % window width in seconds. if this is empty (cfg_in.window = []) it will not use a sliding window.  

cfg = ProcessConfig2(cfg_def, cfg_in); 

%% Filter the signals
if LFP_1.cfg.hdr{1}.SamplingFrequency ~= LFP_2.cfg.hdr{1}.SamplingFrequency
    error('Sampling frequencies between LFP_1 (%dhz) and LFP_2(%dhz) disagree.', LFP_1.cfg.hdr{1}.SamplingFrequency, LFP_2.cfg.hdr{1}.SamplingFrequency)
end

Fs = LFP_1.cfg.hdr{1}.SamplingFrequency; 
cfg_filt = [];
cfg_filt.f = [300 600];
cfg_filt.type = 'butter'; %
cfg_filt.order = 4; % filter order
cfg_filt.display_filter = 0;

% filter LFP 1
Filt_csc_1 = FilterLFP(cfg_filt, LFP_1);

Filt_csc_1 = Filt_csc_1.data; % save some memory by removing extra fields. 


% filter LFP 2
Filt_csc_2 = FilterLFP(cfg_filt, LFP_2);
Filt_csc_2 = Filt_csc_2.data; % save some memory by removing extra fields. 


%% get the xcorr between the two signals using a sliding window. 

% make a sliding window
t_window = round(cfg.window*Fs); 
t_steps = 1:diff(t_window):length(Filt_csc_1); 
csc_corr = zeros(1,length(Filt_csc_1)); 


for iTs = 1:length(t_steps)-1
        csc_corr(t_steps(iTs):t_steps(iTs+1)) = corr(Filt_csc_1(t_steps(iTs):t_steps(iTs+1))', Filt_csc_2(t_steps(iTs):t_steps(iTs+1))'); % get xcorr for this window
end



% put the output into TSD format. 
corr_out = LFP_1;
corr_out.label = 'Corr_EMG';
corr_out.data = csc_corr.*max([LFP_1.data, LFP_2.data]); % normalize to the LFP range




