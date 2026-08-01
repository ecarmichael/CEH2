function HF_TFC_Cond(cfg_in, evts, iv_in); 
%% HF_TFC_Cond:  collects the occurance of events within a head fixed TFC 
% conditioning (typically day 3) for movement, behaviour, and any IV
% inputs.
%
% Inputs: cfg_in [struct]  contains configuration paramters. Any specified
%         will overwrite defaults. 
%
%         evts [struct]  OE events file.  Should contain all of the
%         movement, lick, tones, puffs, ....
%
%         iv_in [struct]  contains the start and stop times of some
%         interval data (typically SWR times]  (optional)

% plot the events as a check. 
figure;
%% initialize

if nargin < 3
    iv_in = []; 
end

cfg_def.start = '5'; 
cfg_def.tone1 = '13'; 
cfg_def.tone2 = '13'; 

cfg_def.mov = '8'; 
cfg_def.lick = '15'; 


%% get the behaviour states