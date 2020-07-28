
data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\Inter';
Subjects = {'PV1043', 'PV1060','PV1069'};


for iSub = 1:length(Subjects)
    cd([data_dir filesep Subjects{iSub}])
    sessions = MS_list_dir_names_any(cd, 'LTD'); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    
    for iSess = 1:length(sessions)
        
        cd([data_dir filesep Subjects{iSub} filesep sessions{iSess}])
        load('all_binary_pre.mat');
        tracker{iSub, iSess}.pre = sum(all_binary_pre, 'all');
        clear all_binary_pre
        
        MS_re_binarize_JC(2, cd, cd, 'ms_resize', 'ms_resize');
        load('all_binary_pre.mat');
        
        tracker{iSub, iSess}.post = sum(all_binary_pre, 'all');
        clear all_binary_pre
    end
end

