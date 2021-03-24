%% sandbox_EV_SWD_SWR_examples

data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\final_analysis'; 
Subs = MS_list_dir_names(data_dir);


Fs = 32000; 


%% collect all the ripples and SWD LFPs for each animal

for iSub = 1:length(Subs)
    sessions = MS_list_dir_names([data_dir filesep Subs{iSub}]);
    
    for iSess = 1:length(sessions)
        
        % get the Fs from a csc file. 
        
        cd([data_dir filesep Subs{iSub} filesep sessions{iSess}])
        load('ripple_1.mat'); 
        ripple_1.LFP{iSub, iSess} = ripple_lfp; 
        ripple_1.data{iSub, iSess} = ripple_data; 
        ripple_1.ind{iSub, iSess} = ripple_ind; 
        
        clear ripple_lfp ripple_data ripple_ind
        
        % get the ripple 2
        load('ripple_2.mat'); 
        ripple_2.LFP{iSub, iSess} = ripple_lfp; 
        ripple_2.data{iSub, iSess} = ripple_data; 
        ripple_2.ind{iSub, iSess} = ripple_ind; 
        
        clear ripple_lfp ripple_data ripple_ind
        
        % get SWD 1
        load('swd_1.mat'); 
        swd_1.LFP{iSub, iSess} = swd_lfp; 
        swd_1.data{iSub, iSess} = swd_data; 
        swd_1.ind{iSub, iSess} = swd_ind; 
        
        clear swd_lfp swd_data swd_ind
        
        % get the ripple 2
        load('swd_2.mat'); 
        swd_2.LFP{iSub, iSess} = swd_lfp; 
        swd_2.data{iSub, iSess} = swd_data; 
        swd_2.ind{iSub, iSess} = swd_ind; 
        
        clear swd_lfp swd_data swd_ind

    end
end
    