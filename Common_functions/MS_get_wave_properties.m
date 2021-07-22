function [wave_prop] = MS_get_wave_properties(S, plot_flag)
%% MS_get_wave_properties:
%
%
%
%    Inputs: 
%    - S [struct]   Spike 'S' structure from LoadSpikes
%
%    - plot_flag: [logical]  0 = no plots, 1 = plots. 
%
%    Outputs: 
%    - wave_prop  [struct] containing:
%                    - burting_idx
%                    - spike_width
%                    - peak_val
%                    - trough_val
%                    - peak2trough
%
%
%
%
% EC 2021-07-22   initial version 
%
%
%
%% initialize