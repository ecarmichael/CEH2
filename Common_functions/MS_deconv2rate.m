function ms_out = MS_deconv2rate(cfg_in, ms_in)
%% MS_deconv2rate: converts devonvoled signals from OASIS into a firing rate following Grosmark et al. 2021.
%
%
%
%    Inputs:
%    - cfg [struct]   configuration see the defaults below.
%
%    - ms_in: [struct] the miniscope struct from cnmfe. If RawTraces is an
%    array it will convert as if the signal is continuous, if it is in
%    cells it will convert each segment/cell.
%
%
%
%    Outputs:
%    - ms_out: [struct]  ms_in but with new 'deconvrate' field.
%
%
%
%
% EC 2022-04-12   initial version
%
%
%
%% initialize
cfg_def = [];

cfg_def.min_decon = 0.01; % minimum to be considered a spike.
cfg_def.bins = 0.05; % bin size for 'spikes'l
cfg_def.gau_win = 1; % window size for gaussian conv
cfg_def.gau_sd = .05; % standard devision for gaussian conv.

cfg = ProcessConfig(cfg_def, cfg_in);
%%

ms_out = ms_in;

if ~iscell(ms_in.deconv)
    fprintf('<strong>%s</strong>: single continuous data detected \n\n', mfilename);
    
    Csp = ms_in.deconv./ms_in.denoise;
    Csp = Csp > cfg.min_decon;
    
    gauss_window = cfg_def.gau_win./cfg_def.bins; % 1 second window
    gauss_SD = cfg_def.gau_sd./cfg_def.bins; % 0.02 seconds (20ms) SD
    gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg_def.bins; % normalize by binsize
    gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window
    
    ms_out.rate = gau_sdf;
    ms_out.rate_cfg = cfg;
    
else
    fprintf('<strong>%s</strong>: segmented data detected (%.0f segments in RawTraces), treating as segmented\n', mfilename, length(ms_in.RawTraces));
    for iSeg = 1:length(ms_in.deconv)
        
        Csp = ms_in.deconv{iSeg}./ms_in.denoise{iSeg};
        Csp = Csp > cfg.min_decon;
        
        gauss_window = cfg_def.gau_win./cfg_def.bins; % 1 second window
        gauss_SD = cfg_def.gau_sd./cfg_def.bins; % 0.02 seconds (20ms) SD
        gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg_def.bins; % normalize by binsize
        gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window
        
        ms_out.rate{iSeg} = gau_sdf;
        
        
    end
    ms_out.rate_cfg = cfg;
    
    
    
end
