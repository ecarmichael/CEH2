%% Subiculum Initialization
%   Run this file to initialize code and data directories. 
%
%
%
%
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;

if ismac
 error('Not setup for macOS yet'); 
elseif strcmp(os, 'GLNXA64')
    
    PARAMS.data_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/'; % where to find the raw data
    PARAMS.inter_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/Inter/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/Stats/'; % where to put the statistical output .txt
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.OASIS_dir = '/home/ecarmichael/Documents/GitHub/OASIS_matlab'; % contains the OASIS decon method: Friedrich et al. 2017: https://doi.org/10.1371/journal.pcbi.1005423
    %
else
    PARAMS.data_dir = 'J:\Williams_Lab\II_classification'; % raw data location.
    PARAMS.inter_dir = 'J:\Williams_Lab\II_classification\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\II_classification\Inter\Stats'; % where to put the statistical output .txt
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.OASIS_dir = '/home/ecarmichael/Documents/GitHub/OASIS_matlab'; % contains the OASIS decon method: Friedrich et al. 2017: https://doi.org/10.1371/journal.pcbi.1005423

end

% Pretty colours
PARAMS.L_grey = [0.8 0.8 0.8];
PARAMS.D_grey = [0.2 0.2 0.2];
PARAMS.blue = [0.3639    0.5755    0.7484];
PARAMS.red = [0.9153    0.2816    0.2878];
PARAMS.green= [0.4416    0.7490    0.4322];
PARAMS.gold = [1.0000    0.5984    0.2000];

rng(11,'twister') % for reproducibility

% add the required code
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

% make the Inter and Stats dir if they don't exist.
if ~exist(PARAMS.stats_dir,'dir')
    mkdir(PARAMS.stats_dir)
end

if ~exist(PARAMS.inter_dir,'dir')
    mkdir(PARAMS.inter_dir)
end

clear d os
beep off % I know when I mess up without that annoying beep, thanks.

%% list all directories for each recording as they are found within the main data directory

Data_dir = {'ck2cre-1359hd/2021_01_30/14_18_06'; ...
    'ck2cre-1359hd/2021_01_30/14_18_06'}; 



