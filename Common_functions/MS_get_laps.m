function [laps, start_idx, end_idx] = MS_get_laps(direction, length_thresh)
%% MS_extract_laps: attempts to pull out linear track laps based on sustained running using a direction vector. 
%
%
%
%    Inputs: 
%    - direction [1 x nSamples]  binary based (1 = running in one
%    direction) 
%
%    - length_thresh [double]  threshold (in samples) to qualify as a lap. 
%
%    Outputs: 
%    - laps [1 x nSamples]  vector with lap number.  (eg:
%    111111110000000022222200000003333333300000...)
%
%    - start_idx [nLaps]   the indicies that correspond to the start of
%    each lap. 
%
%    - end_idx [nLaps]   the indicies that correspond to the end of
%    each lap. 
%
%
% EC 2020-07-15   initial version 
%
%
%
%% initialize

if nargin <2
    error('Speficy a direction vector and a length_threshold (in samples)')
end

%% get the transitions in direction and put them into blocks. 

[pks, loc, width] = findpeaks(direction*1, 'MinPeakWidth', length_thresh);

laps = zeros(size(direction)); % clone size of direction

for iP = 1:length(pks)
    laps(loc(iP):loc(iP)+width(iP)) = iP;
end

start_idx = loc;
end_idx = loc+width; 


