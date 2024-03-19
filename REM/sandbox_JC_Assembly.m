%% sandbox_JC_Assembly

% collect the raw and detrend files;
inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter';
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
            if strcmpi(sess_list(iS).name(end), '3')
                continue
            end
            
            cd([data_dir filesep sub_list(iSub).name filesep sess_list(iS).name]);
            
            
            warning off
            %             load('ms.mat', 'ms')
            %             load('behav.mat')
            %             load('all_binary_pre_REM.mat')
            %             load('all_binary_post_REM.mat')
            
            
            %             info = [];
            %             info.subject = sub_list(iSub).name;
            %             info.session = sess_list(iS).name;
            %             info.nCells = ms.numNeurons;
            %             info.dur_preREM = size(all_binary_pre_REM,1)/30/60;%min
            %             info.dur_postREM = size(all_binary_post_REM,1)/30/60;%min
            
            
            this_file = strsplit(sess_list(iS).name, 'pv');
            
            load([inter_dir filesep 'pv' this_file{2} '_data.mat'])
            
            
            %append the deconv REM data
            
            load('all_detrendRaw_pre_REM.mat')
            load('all_detrendRaw_post_REM.mat')
            load('all_seg_idx.mat');
            
            all_denoise_pre_REM = []; all_deconv_pre_REM = [];
            parfor iChan = 1:size(all_detrendRaw_pre_REM,2)
                tic;
                [denoise,deconv] = deconvolveCa(all_detrendRaw_pre_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
                toc;
                all_denoise_pre_REM(:,iChan) = denoise;    all_deconv_pre_REM(:,iChan) = deconv;
            end
            
            all_denoise_post_REM = []; all_deconv_post_REM = [];
            parfor iChan = 1:size(all_detrendRaw_post_REM,2)
                tic;
                [denoise,deconv] = deconvolveCa(all_detrendRaw_post_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
                toc;
                all_denoise_post_REM(:,iChan) = denoise;    all_deconv_post_REM(:,iChan) = deconv;
            end
            
            
            fprintf('Saving %s   %s ...\n', info.subject, info.session)
            save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav',...
                'all_deconv_pre_REM', 'all_deconv_post_REM',...
                'all_denoise_pre_REM', 'all_denoise_post_REM',...
                'all_binary_pre_REM', 'all_binary_post_REM',...
                'all_seg_idx', 'info')
            warning on
        end
    end
end



%%

cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')


load('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\pv1043_LTD1_data.mat')

% MS_PCA_ICA_no_fig(

% MS_PCA_ICA_CA([], ms, behav)

%%

inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Inter';

cd(inter_dir)
f_list = dir('*data.mat');

for iF = length(f_list):-1:1
    
    load([inter_dir filesep f_list(iF).name]);
    
    sub = info.subject;
    sess = info.session;
    
    cd([data_dir filesep lower(sub)])
    
    s_dir = dir(['*' sess]);
    
    cd([data_dir filesep upper(sub) filesep s_dir.name]);
    
    
    %append the deconv REM data
    
    load('all_detrendRaw_pre_REM.mat')
    load('all_detrendRaw_post_REM.mat')
    load('all_seg_idx.mat');
    
    all_denoise_pre_REM = []; all_deconv_pre_REM = [];
    parfor iChan = 1:size(all_detrendRaw_pre_REM,2)
        tic;
        [denoise,deconv] = deconvolveCa(all_detrendRaw_pre_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
        toc;
        all_denoise_pre_REM(:,iChan) = denoise;    all_deconv_pre_REM(:,iChan) = deconv;
    end
    
    all_denoise_post_REM = []; all_deconv_post_REM = [];
    parfor iChan = 1:size(all_detrendRaw_post_REM,2)
        tic;
        [denoise,deconv] = deconvolveCa(all_detrendRaw_post_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
        toc;
        all_denoise_post_REM(:,iChan) = denoise;    all_deconv_post_REM(:,iChan) = deconv;
    end
    
    
    fprintf('Saving %s   %s ...\n', info.subject, info.session)
    save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav',...
        'all_deconv_pre_REM', 'all_deconv_post_REM',...
        'all_denoise_pre_REM', 'all_denoise_post_REM',...
        'all_binary_pre_REM', 'all_binary_post_REM',...
        'all_seg_idx', 'info')
    warning on
    
    clear all_* ms behav info
end

%%
inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Inter';

cd(inter_dir)
f_list = dir('*data.mat');

for iF = length(f_list):-1:1
    
    load([inter_dir filesep f_list(iF).name]);
    
    sub = info.subject;
    sess = info.session;
    
    cd([data_dir filesep lower(sub)])
    
    s_dir = dir(['*' sess]);
    
    cd([data_dir filesep upper(sub) filesep s_dir.name]);
    
    
    %append the deconv REM data
    
    load('all_detrendRaw_pre_REM.mat')
    load('all_detrendRaw_post_REM.mat')
   
    
    fprintf('Saving %s   %s ...\n', info.subject, info.session)
    save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav',...
        'all_deconv_pre_REM', 'all_deconv_post_REM',...
        'all_denoise_pre_REM', 'all_denoise_post_REM',...
        'all_binary_pre_REM', 'all_binary_post_REM',...
        'all_detrendRaw_pre_REM', 'all_detrendRaw_post_REM',...
        'all_seg_idx', 'info')
    warning on
    
    clear all_* ms behav info
end