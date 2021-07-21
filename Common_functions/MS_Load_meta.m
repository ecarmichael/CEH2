function Meta = MS_Load_meta(file_name)
%% MS_Load_meta: loading function for the '*meta.m' files.  Meta files are 
% used to track all experimental variables like subject, tragets, behaviour...
% 
%
%
%   Inputs:
%       - none: will look for anything ending in '*Meta.m'
%       - file_name: [string] file to run as a string (eg:
%       'M21_2021_07_06_OF_meta.m', or
%       '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inProcess/M21/M21_2021_07_06_OF_meta.m').
%        
%   Outputs:
%       - none.  Instead it will run the *meta.m file resulting in the meta
%       variable. 
%
%
%   EC 2019-12-30: initial version based on ExpKeys/LoadExpKeys/FindFile by
%   MvdM and AD Redish
%
%   EC 2021-07-16 updated to work for new 'meta.m' files to avoid stepping
%   on Redish/vdm lab style. 
%
%% extract inputs and run or 

if nargin ==1 % if specified, just run it. 
    run(file_name); 
    fprintf(['\n<strong>%s</strong>: ' file_name ' loaded.\n'], mfilename); 
    
else
    file_name = dir('*meta.m');
    if isempty(file_name)
        error('No *meta.m files found here')
    end
        run(file_name.name); 
    fprintf(['\n<strong>%s</strong>: ' file_name.name ' loaded.\n'], mfilename); 
    
end
clear file_name
end