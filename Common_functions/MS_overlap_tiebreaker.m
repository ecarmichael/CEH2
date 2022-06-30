function [idx_out] = MS_overlap_tiebreaker(data_in, val_in, win_s)
%% MS_tiebreaker: identifies overlapping values (within n samples from earch other) in data_in and finds the value with the highest corresponding val_in value
%   Example:  data_in = [1   5    6   7   12]
%             val_in  = [.2  .4   .6  .1  .4]
%             win_s = 3
%              .... 5, 6, and 7 are within the overlap, the tie break goes
%              to 6 
%             idx_out = [1, 6, 12]
%
%    Inputs: 
%    - data_in: [nSamples x 1]    data to check for overlaps (<= win_s)
%
%    - val _in: [nSamples x 1]    values for breaking the tie. must be
%    equal size to data_in
%
%    - win_s: [double]                number of samples for the overlap window
%
%    Outputs: 
%    -
%
%
%
%
% EC 2022-06-03   initial version 
%
%
%
%% intialize

if nargin < 3
    warning('not enough inputs. Running demo...')
    data_in = [1   5    6   7   12];
    val_in  = [.2  .4   .6  .1  .4];
    win_s = 3 ;
end

if length(data_in) ~= length(val_in)
    error('data_in and val_in differ in size')
end

%% convert data to start and ends
data = [data_in, data_in+win_s];









 