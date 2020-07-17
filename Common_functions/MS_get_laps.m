function [laps, start_idx, end_idx] = MS_get_laps(direction, length_thresh, merge_width)
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
%    - merge_width [double] 
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
elseif nargin ==2 
    merge_width = 0; 
end

%% get the transitions in direction and put them into blocks. 

[pks, loc, width] = findpeaks(direction*1, 'MinPeakWidth', length_thresh);

laps = zeros(size(direction)); % clone size of direction

% get the diff between peaks in samples (used for mergers)
peak_diff = diff(loc);
merge_idx = peak_diff < merge_width; 

% get the start and end of each block. 
start_idx = loc;
end_idx = loc+width; 

for iP = 1:length(pks)
    laps(loc(iP):loc(iP)+width(iP)) = iP;
end

% merge peaks within a certain distance
if merge_width ~=0
    laps_n = 2;
    for iM =2:length(pks)
        if merge_idx(iM-1)
            laps(loc(iM-1):loc(iM)+width(iM)) = laps(loc(iM-1));
            start_idx(iM-1) = loc(iM-1);
            end_idx(iM-1) = loc(iM)+width(iM);
            start_idx(iM) = NaN;
            end_idx(iM) = NaN;
        else
            laps(loc(iM):loc(iM)+width(iM)) = laps_n;
            laps_n = laps_n+1;
        end
    end
    
% remove the merged start 
    start_idx(isnan(start_idx)) = [];
    end_idx(isnan(end_idx)) = [];
end




