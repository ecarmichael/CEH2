%% CE_radial_preprocess


sub_list = {'1799', '1842', '1860', 'dSub_g8_E', 'dSub_g8_F'};

if ispc
    usr_path = ['C:\Users\' getenv('username')];
end

data_dir = [ usr_path  strrep('\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Rad_maze', '\', filesep)];
inter_dir = [data_dir filesep 'Inter'];

cd(data_dir)


%% get all the raw data


for iSub = 1:length(sub_list)
    
    
    s_list = dir([data_dir filesep sub_list{iSub} filesep '20*']);
    
    % loop over sessions
    for iS = 1:length(s_list)
        ms_check =[]; keep_check = []; enc_check = []; rec_check = [];
        
        fprintf([sub_list{iSub} '-' s_list(iS).name ':  '])
        if exist([inter_dir filesep sub_list{iSub} filesep s_list(iS).name filesep 'ms.mat'], 'file')
            ms_check = 1;
        else
            ms_check = 0;
        end
        fprintf('ms.mat  %d   |   ', ms_check)
        
        % check for keep_index
        if exist([inter_dir filesep sub_list{iSub} filesep s_list(iS).name filesep 'keep_idx.mat'], 'file')
            keep_check = 1;
        else
            keep_check = 0;
        end
        fprintf('keep_idx.mat  %d   |   ', keep_check)
        
        
        % check for enc behav
        if exist([inter_dir filesep sub_list{iSub} filesep s_list(iS).name filesep 'behav_enc.mat'], 'file')
            enc_check = 1;
        else
            enc_check = 0;
        end
        fprintf('behav_enc.mat  %d   |   ', enc_check)
        
        if exist([inter_dir filesep sub_list{iSub} filesep s_list(iS).name filesep 'behav_enc.mat'], 'file')
            rec_check = 1;
        else
            rec_check = 0;
        end
        fprintf('behav_rec.mat  %d   |   ', rec_check)
        
        if sum([ms_check keep_check enc_check rec_check]) == 4
            fprintf('done')
        end
        
        fprintf('\n')
        
    end
    
    
end