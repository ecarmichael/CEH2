function mua = MS_MUA(S, tvec, sigma)
%% MS_MUA:
%
%
%
%    Inputs: 
%    - S [struct]  Spikes in the TS format
%    - tvec [1xnsamples]  time vector
%
%
%    Outputs: 
%    - mua [struct]  
%
%
%
%
% EC 2026-03-06   initial version 
%  based on rawMUA by Aacarey https://github.com/vandermeerlab/vandermeerlab/blob/d23da821fe75f87c5e54d62b1c2e1f9aa37cfa60/code-matlab/tasks/Alyssa_Tmaze/beta/rawMUA.m#L49
%
%
%% initialize

if nargin < 3
    sigma = 1; % Standard deviation for Gaussian kernel
end

% combine spikes
all_S = []; 
for ii = 1:length(S.t)
 all_S = [all_S; S.t{ii}]; 

end

% sort in time. 
all_S = sort(all_S); 

all_S_bin = zeros(size(tvec)); 
for ii = 1:length(all_S)
    this_s = nearest_idx3(all_S(ii), tvec);
    all_S_bin(this_s) = all_S_bin(this_s)+1;

end

% convert to time series data
mua  = tsd(tvec, all_S_bin'); 
% smoothing

smoothed_MUA = imgaussfilt(mua.data, sigma);
mua.data = smoothed_MUA;