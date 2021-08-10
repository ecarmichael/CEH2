function [f_names, f_fold] = MS_list_dir_names(dir_name, str_to_find)
%% MS_list_dir_names: get the names of all the files/folders in the dir_name. Removes the '.' and '..' dirs and can look for a specific string 'str_to_find'
%    MS_list_dir_names returns the dir name and dir of all dirs that
%    contain all of the str_to_find.  For cases where any of the str_to_find are
%    found use MS_list_dir_names_any
%
%
%    Inputs:
%     - dir_name: [string] path to dir to list.  If empty will use current
%     dir.
%
%     - str_to_find: [n x 1 cell array] strings to look for in the dir.
%           example: find dir with {'day1', 'base1'}
%                   [ 0   0   1   1   0   0   0   0]
%                   [ 0   0   1   0   1   0   1   0]
%                   returns
%      keep_idx =   [ 0   0   1   0   0   0   0   0]
%
%    Outputs:
%     - f_name: [n x1 cell array] contains the names of the folders/files
%     in the dir.
%     - f_full: [n X 1 cell array] contains the full folder path.
%
%
% EC 2020-03-02   initial version
%
% See also MS_list_dir_names

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
this_dir = dir(dir_name);

% switch between using the str_to_find
if  ~isempty(str_to_find)
    %workaround to get str_to_find into a cell array. 
    if ischar(str_to_find) == 1
        str_to_find = {str_to_find};
    end
    % cycle through strings to find.
    for iStr = 1:length(str_to_find)
        for iF = 1:length(this_dir)
            if strcmpi(this_dir(iF).name, '.') || strcmpi(this_dir(iF).name, '..')
                continue
            elseif contains(lower(this_dir(iF).name), lower(str_to_find{iStr}))
                this_f_names{iStr}{iF} = this_dir(iF).name;
                this_f_fold{iStr}{iF} = this_dir(iF).folder;
            else
                this_f_names{iStr}{iF} = [];
                this_f_fold{iStr}{iF} = [];
            end
        end
        % keep the cases where all of the required strings were found.
        keep_mat(iStr,:) = ~cellfun(@isempty,this_f_names{iStr});
    end
    
    % get all the cases where all the str_to_find are all present in the
    % dir name.
    f_names_found = {this_dir(sum(keep_mat,1) == length(str_to_find)).name};
    f_fold_found = {this_dir(sum(keep_mat,1) == length(str_to_find)).folder};
else
    for iF = 1:length(this_dir)
        if strcmp(this_dir(iF).name, '.') || strcmp(this_dir(iF).name, '..')
            continue
        else
            f_names_found{iF} = this_dir(iF).name;
            f_fold_found{iF} = this_dir(iF).folder;
        end
    end
    f_names_found = f_names_found(~cellfun('isempty',f_names_found));
    f_fold_found = f_fold_found(~cellfun('isempty',f_fold_found));
end

%% quick clean up

% if only one dir found, convert it to a cell array for output. 
if ischar(f_names_found)
    f_names{1} = f_names_found;
    f_fold{1} = f_fold_found;
else
    f_names = f_names_found;
    f_fold = f_fold_found; 
end
