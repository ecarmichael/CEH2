function [pos, behav] = MS_DLC2TSD_single(fname, conv_fac, plot_flag)
%% MS_DLC2TSD: loads and collects all DLC files in a directory. Will skip over files without a number since DLC saves the interation number in the .csv
%
%
%    Inputs:
%    -  dir_in:  [string]  name of the directory to collect.
%
%    - plot_flag: bool do you want to see a movie of the HD / position?
%
%    Outputs:
%    - pos [struct] position data in the tsd format.
%
%    - behav [struct]  data in the 'behav' format for miniscope analysis.
%    [defaults is 0]
%
%
%
%
% EC 2023-01-28  initial version
%% initialize

if nargin == 1
    error('needs a video file for time');
elseif nargin == 2
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==3
    plot_flag = 0;
end