%% sandbox_convert_EV_data

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\RAW Calcium\Inter' ;
behav_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\final_analysis';
inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV';


sub_list = {'535', '537', '540'};

for iS = 1:length(sub_list)
    cd(data_dir);
    
    d = dir(['*' sub_list{iS} '*']);
    
    
    cd([data_dir filesep  d.name])
    
    d_list = dir('*day*');
    
    for iD = 1:length(d_list)
        cd([d_list(iD).folder filesep  d_list(iD).name])
        
        load('ms_trk.mat');
        ms = ms_trk; clear ms_trk;
        
        load('all_binary_pre_REM.mat')
        load('all_binary_post_REM.mat')
        
        load('all_seg_idx.mat');
        
        info = [];
        info.subject = sub_list{iS};
        info.session = d_list(iD).name;
        info.nCells = ms.numNeurons;
        info.dur_preREM = size(all_binary_pre_REM,1)/30/60;%min
        info.dur_postREM = size(all_binary_post_REM,1)/30/60;%min
        
        % move over to the behaviour folder.
        cd([behav_dir filesep sub_list{iS} filesep d_list(iD).name])
        
        load('behav.mat');

        fprintf('Saving %s   %s ...\n', info.subject, info.session)
        save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav', 'all_binary_pre_REM', 'all_binary_post_REM','all_seg_idx', 'info')
        
        
    end
end