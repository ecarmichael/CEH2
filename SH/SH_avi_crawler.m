function SH_avi_crawler(data_dir)
%% SH_avi_conv_crawler: looks for all directories within a parent for .avi files begining with a specific prefix and runs the ffmpeg conversion to get the 'grey' avi.
%
%    Runs the following shell script which needs to be saved as an alias in
%    the ~/.bashrc.
%
%
%
%
%
%    Inputs:
%    - data_dir [path]   parent directory to be searched.
%
%    - prefix [string]  pattern to look for. Default is any numeric value since this is
%    how the new UCLA miniscope software saves .avi files. Old software
%    uses 'ms'.  NOT USED ATM. 
%
%    Outputs:
%    - none
%
%
%
%
% EC 2023-06-19   initial version
%
%
%
%% initialize
avi_conv_fname = which('avi_conv.sh');

if ~isempty(avi_conv_fname)
    fprintf('<strong>%s</strong>: using avi_conv.sh from <strong>%s</strong>: \n', mfilename, avi_conv_fname);
else
    error('Could not find avi_conv.sh');
end


og_dir = pwd;
dir_list = dir(data_dir);
dir_list(1:2) = []; % remove the first two dirs which are '.' and '..';


fprintf('<strong>%s</strong>: directories to be processed: \n', mfilename);

%% conver to avi grey using avi_conv

for ii = 1:length(dir_list)
    fprintf('%s%s<strong>%s</strong> \n',dir_list(ii).folder, filesep,  dir_list(ii).name)
    
    Scope_dir = dir (fullfile([dir_list(ii).folder filesep dir_list(ii).name], '**', '*Mini*'));
    cd([dir_list(ii).folder filesep dir_list(ii).name filesep Scope_dir(1).name])
    
    % run a shell script for avi conversion using ffmpeg. 
    cmstr = sprintf('%s', ['sh ' avi_conv_fname]);
    system(cmstr);
    
    
end

fprintf('<strong>%s</strong>: all avi files processed: \n', mfilename);

% bounce back to original dir
cd(og_dir); 
