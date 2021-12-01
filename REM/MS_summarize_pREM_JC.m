function MS_summarize_pREM_JC(data_dir)
%% MS_summarize_pREM_JC:
%
%
%
%    Inputs: 
%    - data_dir [path]   path to processed pREM data.  Script will group
%          data based on file structure.  Eg if 'pREM_dur.mat' is found in 
%          'PV1043/6_11_2019_PV1043_LTD1' then it will be grouped as
%          PV1043. 
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2021-09-22   initial version 

%% initliaze

if nargin < 1
    data_dir = cd; 
end

%% look for pREM_dur files in data_dir.

pREM_files = dir(fullfile(data_dir,'**', 'pREM_dur.mat')); % find all instances of the pREM_dur.mat files
for iF = length(pREM_files):-1:1
   pREM_paths{iF} = [pREM_files(iF).folder filesep pREM_files(iF).name]; % convert dir output to just full paths with filenames. 
   pREM_paths{iF} = pREM_paths{iF}(end-strfind(reverse(pREM_paths{iF}), reverse(data_dir))+3:end); % get only the dir paths beyond the data_dir; 
end

fprintf('<strong>%0.0f</strong> pREM files detected in <strong>%s</strong>\n', length(pREM_paths), data_dir)

%% loop over pREM files found and save the values in groups if needed
subjects = dir(data_dir); subjects  = {subjects(3:end).name}; % get a list of subjects

for iP = length(pREM_paths):-1:1
    parts = strsplit(pREM_paths{iP}, filesep);
    load([data_dir filesep pREM_paths{iP}], 'pREM_dur'); 
    
    mean_out.(parts{1}).(['D_' parts{2}]).REM_all = pREM_dur.mean_REM_prct_all;
    mean_out.(parts{1}).(['D_' parts{2}]).REM_pre = pREM_dur.mean_REM_prct_pre; 
    mean_out.(parts{1}).(['D_' parts{2}]).REM_post = pREM_dur.mean_REM_prct_post; 

end






