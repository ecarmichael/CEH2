%% sandbox_convert_EV_data

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\RAW Calcium\Inter' ;
behav_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\final_analysis';
inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV';
oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';

cd(oasis_dir)
addpath(genpath(oasis_dir));
oasis_setup

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
        
        ms.time = ms.time - ms.time(1); 
       
        ms.RawTraces = ms.RawTraces - mean(ms.RawTraces);
        
        if ~isfield(ms, 'deconv') % get the deconvolved trace if not already present.
            ms = MS_append_deconv(ms, 1);
        end
            
        load('all_binary_pre_REM.mat')
        load('all_binary_post_REM.mat')
        
                
        load('all_detrendRaw_pre_REM.mat')
        load('all_detrendRaw_post_REM.mat')
        
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
        
        load('all_seg_idx.mat');
        
        info = [];
        info.subject = sub_list{iS};
        info.session = ['LTD' d_list(iD).name(end)];
        
        
        info.nCells = ms.numNeurons;
        info.dur_preREM = size(all_binary_pre_REM,1)/30/60;%min
        info.dur_postREM = size(all_binary_post_REM,1)/30/60;%min
        
        % move over to the behaviour folder.
        cd([behav_dir filesep sub_list{iS} filesep d_list(iD).name])
        
        load('behav.mat');
        
        nan_idx = (behav.position(:,1) < 0 ) | (behav.position(:,1) > 100);
        behav.position(nan_idx,1) = NaN; 
        behav.position(:,1) = fillmissing(behav.position(:,1), 'nearest');
        
        nan_idx = (behav.position(:,2) < 0 );
        behav.position(nan_idx,2) = NaN; 
        behav.position(:,2) = fillmissing(behav.position(:,2), 'nearest');

        fprintf('Saving %s   %s ...\n', info.subject, info.session)
        save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav',...
            'all_binary_pre_REM', 'all_binary_post_REM',...
            'all_deconv_pre_REM', 'all_deconv_post_REM',...
            'all_denoise_pre_REM', 'all_denoise_post_REM',...
            'all_seg_idx', 'info')
        
        
    end
end

%% append deconv
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
cd(oasis_dir)
addpath(genpath(oasis_dir));
oasis_setup

cd(inter_dir)
f_list = dir('*data.mat');

for ii = 1:length(f_list)
    
    load(f_list(ii).name)
    
    fprintf('Processing %s   %s ...\n', info.subject, info.session)
    
    ms.RawTraces = ms.RawTraces - mean(ms.RawTraces); 
    
    if ~isfield(ms, 'deconv') % get the deconvolved trace if not already present.
    ms = MS_append_deconv(ms, 1);
    end
    
    fprintf('Saving %s   %s ...\n', info.subject, info.session)
    save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav', 'all_binary_pre_REM', 'all_binary_post_REM','all_seg_idx', 'info')
    
    
end

%%
data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eva\RAW Calcium\Inter' ;
inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV';


cd(inter_dir); 

f_list = dir('*data.mat');

for iS = 1:length(f_list)

    load([inter_dir filesep f_list(iS).name])

    this_dir = dir([data_dir filesep info.subject filesep '*' info.session(end)]);
    
    load([this_dir.folder filesep this_dir.name filesep 'all_detrendRaw_pre_REM.mat'])
    load([this_dir.folder filesep this_dir.name filesep 'all_detrendRaw_post_REM.mat'])

        
        fprintf('Saving %s   %s ...\n', info.subject, info.session)
        save([inter_dir filesep info.subject '_' info.session '_data.mat'], 'ms', 'behav',...
            'all_binary_pre_REM', 'all_binary_post_REM',...
            'all_deconv_pre_REM', 'all_deconv_post_REM',...
            'all_denoise_pre_REM', 'all_denoise_post_REM',...
            'all_detrendRaw_pre_REM', 'all_detrendRaw_post_REM',...
            'all_seg_idx', 'info')
        
        clear all_* ms behav info
end
