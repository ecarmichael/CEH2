function [wave_prop] = MS_get_wave_properties(S, wave,tvec,  plot_flag)
%% MS_get_wave_properties:
%
%
%
%    Inputs: 
%    - S [struct]   Spike 'S' structure from LoadSpikes
%
%    - wave  [5 x nSamples]  waveform array. 
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
if nargin < 4
    plot_flag = 0;
end

%% 

spike_rate = length(S.t{1})./(tvec(end) - tvec(1)); 


%% get the autocorrelation (from vandermeerlab.org wiki: https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:nsb2016:week9)
xbin_centers = -.025-0.001:0.001:0.025+0.001; % first and last bins are to be deleted later
ac = zeros(size(xbin_centers));
 
for iSpk = 1:length(S.t{1})
 
   relative_spk_t = S.t{1} - S.t{1}(iSpk);
 
   ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.
 
end
 
xbin = xbin_centers(2:end-1); % remove unwanted bins
zero_idx = find(xbin == 0); 
ac = ac(2:end-1);
 ac(zero_idx) = 0; 

% ac = ac./max(ac); % normalize

figure(900)
bar(xbin*1000, ac)
xlabel('Lag (ms)');


bin_idx_2_6 = zero_idx +2: zero_idx+6; 
bin_idx_2_20 = zero_idx +2: zero_idx+20; 

