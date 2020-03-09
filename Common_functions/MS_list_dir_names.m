function [f_names, f_fold] = MS_list_dir_names(dir_name, str_to_find)
%% MS_list_dir_names: get the names of all the files/folders in the dir_name. Removes the '.' and '..' dirs and can look for a specific string 'str_to_find'
%
%
%
%    Inputs:
%     - dir_name: [string] path to dir to list.  If empty will use current
%     dir.
%
%     - str_to_find: [n x 1 cell array] strings to look for in the dir.
%           example:
%
%    Outputs:
%     - f_name: [n x1 cell array] contains the names of the folders/files
%     in the dir.
%     - f_full: [n X 1 cell array] contains the full folder path. 
%
%
%
%
% EC 2020-03-02   initial version
%
%
%% initialize

if nargin == 0
    dir_name = cd;
    str_to_find = [];
    fprintf('<strong>%s</strong>:, no dir_name or str_to_find inputs specified. Using current folder and no str_to_find\n', mfilename);
elseif nargin == 1
    str_to_find = [];
    fprintf('<strong>%s</strong>:,using dir: %s. no str_to_find inputs specified. Usin no str_to_find\n', mfilename, dir_name);
end

%% get the name
f_names = [];
this_dir = dir(dir_name);

% switch between using the str_to_find
if  ~isempty(str_to_find)
    
    for iF = 1:length(this_dir)
        if strcmp(this_dir(iF).name, '.') || strcmp(this_dir(iF).name, '..')
            continue
        elseif find(contains(this_dir(iF).name, str_to_find))
            f_names{iF} = this_dir(iF).name;
            f_fold{iF} = this_dir(iF).folder; 
        else
            f_names{iF} = [];
            f_fold{iF} = []; 
        end
    end
    
else
    for iF = 1:length(this_dir)
        if strcmp(this_dir(iF).name, '.') || strcmp(this_dir(iF).name, '..')
            continue
        else
            f_names{iF} = this_dir(iF).name;
            f_fold{iF} = this_dir(iF).folder; 
        end
    end
end


f_names = f_names(~cellfun('isempty',f_names));
f_fold = f_fold(~cellfun('isempty',f_fold));

