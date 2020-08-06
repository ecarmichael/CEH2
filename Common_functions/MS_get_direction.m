function idx_out = MS_get_direction(position, threshold, method)
%% MS_get_direction: gets the indicies of movement direction from a 1d position vector. [Default]Use -threshold for leftbound and + for right.
%
%
%
%    Inputs:
%    - position [1 x nSamples] position vector.
%
%    - threshold [double]  cutoff threshold for movement distance traveled.
%    by default - values are for leftbound amd + are for right bound.
%    Depends on the units and camera setup. Units based on the method
%
%    - method [string]  can be 'zcore' for zscore values, or 'raw' for raw
%    values.  Raw values are the same as the units in the position vector.
%    [Default] 'raw'.
%
%    Outputs:
%    - idx [1 x nSamples] index of all samples exceeding the threshold.
%
%
%
%
% EC 2020-07-15   initial version based on isolate_direction by Guillame
% Etter: https://github.com/etterguillaume/CaImDecoding/blob/master/isolate_direction.m
%
%
%% initialize

if nargin ==1
    error('no threshold specified')
elseif nargin == 2
    method = 'raw'; % default method is 'raw'
end

%%  get the change in position per sample

diff_pos = diff(position);
diff_pos(end+1) = 0; % correct for diff offset.

if strcmpi(method, 'zscore')
    diff_pos = zscore(diff_pos);
end

if threshold < 0
    
    idx_out = diff_pos < threshold;
    
elseif threshold > 0
    
    idx_out = diff_pos > threshold;
    
else
    error('threshold is zero')
end
end

