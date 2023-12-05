%% sandbox_JC_Assembly

% collect the raw and detrend files;

cd(data_dir)
sub_list = dir('pv*');

for iSub = length(sub_list):-1:1
    cd([data_dir filesep sub_list(iSub).name]);
    
    
    for iSub = length(sub_list):-1:1
        cd([data_dir filesep sub_list(iSub).name]);
        
        
        LT_list = dir('*LTD*');
        HAT_list= dir('*HATD*');
        sess_list=[LT_list; HAT_list];
        
        
        for iS  = length(sess_list):-1:1
            
            cd([data_dir filesep sub_list(iSub).name filesep sess_list(iS).name]);
            
            
            warning off
            load('ms.mat', 'ms')
            load('behav.mat')
            load('all_binary_pre_REM.mat')
            load('all_binary_post_REM.mat')
            
            
            info = [];
            info.subject = sub_list(iSub).name;
            info.session = sess_list(iS).name;
            info.nCells = ms.numNeurons;
            info.dur_preREM = size(all_binary_pre_REM,1)/30/60;%min
            info.dur_postREM = size(all_binary_post_REM,1)/30/60;%min
            
            
            fprintf('Saving %s   %s ...\n', info.subject, info.session)
            save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav', 'all_binary_pre_REM', 'all_binary_post_REM', 'info')
            warning on
        end
    end
end



%%

cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')


load('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\pv1043_LTD1_data.mat')

MS_PCA_ICA_no_fig(

% MS_PCA_ICA_CA([], ms, behav)