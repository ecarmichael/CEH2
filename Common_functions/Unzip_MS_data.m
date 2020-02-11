function Unzip_MS_data(data_dir, target_dir)
%% unzip a file into a folder of the same name at the target_dir.  
%
%
%
%
% EC 2020-02-11   initial version 
%
%
%%

% target_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5';
if isempty(data_dir)
    these_files = dir(fullfile(pwd, '*.zip'));
else
    these_files = dir(fullfile(data_dir, '*.zip'));
end

for iF = 1:length(these_files)
    
    if ~strcmp(these_files(iF).name(1), 'H')
        fprintf('Skipping %s...\n', these_files(iF).name)
        continue
    else
        fprintf('Unziping %s...\n', these_files(iF).name)
        unzip(these_files(iF).name, [target_dir filesep these_files(iF).name(1:end-4)])
        
    end
end