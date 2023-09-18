function [ppx,f] = MS_get_psd(cfg_in, csc)
%% MS_get_psd:
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-09-17   initial version 
%
%
%
%% initialize


cfg_def = [];
cfg_def.hann_win = 2^12;

cfg = ProcessConfig(cfg_def, cfg_in);

%%


[ppx,f] = pwelch(csc.data, hanning(cfg.hann_win), cfg.hann_win/2, cfg.hann_win*2 , csc.cfg.hdr{1}.SamplingFrequency);
