function [laps] = MS_get_laps(direction, length_thresh)
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
%
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

findpeaks(diff(direction))
