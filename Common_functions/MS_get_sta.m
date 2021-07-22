function [st_mat, t_win] = MS_get_sta(cfg_in, csc, S)
%% MS_get_sta: get the spike triggered average between the LFP(csc) and spike times (S.t)
%
%
%
%    Inputs: 
%    - cfg_in [struct] configuration structure.  unused fields will be
%    overwritten by defaults. 
%
%    - csc [struct]  'continuously sampled channel'  LFP data in the 'csc'
%    format from MS_LoadCSC
%
%    - S [struct]  'Spike' stucture containing spike times from
%    LoadSpikes.m
%
%
%
%    Outputs: 
%    - st_mat  [nSamples x nSpikes] spike-triggered LFP matrix to be
%    averaged to get the spike-triggered average. 
%
%    - t_win = time vector corresponding to the st_mat; 
%
%
% EC 2021-07-19   initial version based on vanderMeer lab data analysis
% wiki: https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:nsb2016:week13
%
%
%
%% initialize
cfg_def = [];
cfg_def.t_win = [-1 1]; % time window to use. 
cfg_def.S_chan = 1; % which spike to use. 
cfg_def.csc_chan = 1; % which csc channel to use. 

cfg = ProcessConfig(cfg_def, cfg_in);

%% detrend the data if chronux is present
if exist('locdetrend', 'file')
    
    csc.data(csc.chan,:) = locdetrend(csc.data(cfg.csc_chan,:), csc.cfg.hdr{1}.SamplingFrequency, [1 0.5]);
else
%     [bFilt,aFilt] = butter(2,  2/(csc.cfg.hdr{1}.SamplingFrequency), 'low');
%     fvtool(bFilt, aFilt); 
%     csc.data(csc.chan,:) = filtfilt(bFilt, aFilt, csc.data(cfg.csc_chan,:));
csc.data(cfg.csc_chan,:) = detrend(csc.data(cfg.csc_chan,:));
end

%% generate a time window
t_win = cfg.t_win(1):1/csc.cfg.hdr{1}.SamplingFrequency: cfg.t_win(2);

S_t = S.t{cfg.S_chan};

h = waitbar(0,sprintf('Cell %d/%d...',cfg.S_chan,length(S_t)));

st_mat = [];tic
for iS = length(S_t):-1:1
    
    sta_t = S_t(iS)+t_win(1);
%     if ~isunix
%     NOTE if you do not have a compiled mex file 
        sta_idx = nearest_idx3(csc.tvec,sta_t); % find index of leading window edge
%     else
%         sta_idx = nearest_idx(sta_t, csc.tvec); % find index of leading window edge
%     end
    if sta_idx+length(t_win)-1 > length(csc.data(cfg.csc_chan,:)) % spike spikes where the post-spike window exceeds the csc length
        continue
    else
        this_st_mat = csc.data(sta_idx:sta_idx+length(t_win)-1); % grab LFP snippet for this window
        st_mat(iS,:) = this_st_mat';
    end
   waitbar(iS/length(S_t));

end
toc
close(h) % close the waitbar

