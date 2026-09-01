function S_out = HF_deep_super(S, csc, swr_iv, probe, chan_idx)

%% HF_deep_super: uses the SWR times to create a ripple average for all the channels on a shank by shank basis. Determines the middle of the ripple using 
%
%
%
%    Inputs: 
%    - S [struct]  spike times in the TS format
%
%    - csc [struct]  LFP in the tsd format
%
%    - swr_iv [struct]  swr events in the iv format
%
%    - probe [string]  can be either 'A4x16', A5x12', 'Buz32'
%
%    - chan_idx [n x 1] indices of channels to use. 
%
%    Outputs: 
%    - S_out [struct]  spike TS with the deep/super classification
%
%
%
%
% EC 2026-09-01   initial version 
%
%
%
%% initialize

if nargin < 5
    chan_idx = size(csc.data,1);
end

if strcmpi(probe, 'A4x16')





end


