clear all
close all

ms_dir = {'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_24_2021_pv1252_HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_14_2021_pv1252_LTD3',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_12_2019_PV1069_LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch'};


plot_f = 1; % toggle plots for speed.

cutoff_n_cell = 50;

method = 'decon'; 

for iSess = 4%:length(sessions)
    %% load some data
    
    cd(ms_dir{iSess})
    load('ms_resize.mat');
       load('ms_trk.mat'); 

if strcmpi(method, 'decon')
    load('all_binary_pre_REM.mat');
    load('all_binary_post_REM.mat');
    load('all_binary_pre_SW.mat');
    load('all_binary_post_SW.mat');
    
    data_pre_REM = all_binary_pre_REM; 
    data_post_REM = all_binary_post_REM; 
    
    data_pre_REM = all_binary_pre_REM; 
    data_post_REM = all_binary_post_REM; 
    
    data_trk = ms_trk.Binary; 
    clear all_*

else
    addpath('C:\Users\ecarm\Documents\GitHub\OASIS_matlab')
    oasis_setup;
    
    load('all_RawTraces_pre_REM.mat');
    load('all_RawTraces_post_REM.mat');
    
    load('all_RawTraces_pre_SW.mat');
    load('all_RawTraces_post_SW.mat');
    
    % decon the traces
    all_denoise = []; all_deconv= [];
    fprintf('\n<strong>%s</strong>: deconvolving traces...\n', mfilename)
    tic;
    for iChan = size(all_RawTraces_post_REM,2):-1:1
        disp(num2str(iChan))
        [denoise,deconv] = deconvolveCa(all_RawTraces_pre_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
        ms_pre_REM.denoise(:,iChan) = denoise;    ms_pre_REM.deconv(:,iChan) = deconv;
        
        [denoise,deconv] = deconvolveCa(all_RawTraces_post_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
        ms_post_REMs.denoise(:,iChan) = denoise;    ms_post_REM.deconv(:,iChan) = deconv;
        
        [denoise,deconv] = deconvolveCa(all_RawTraces_pre_SW(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
        ms_pre_SW.denoise(:,iChan) = denoise;    ms_pre_SW.deconv(:,iChan) = deconv;
        
        [denoise,deconv] = deconvolveCa(all_RawTraces_post_SW(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
        ms_post_SW.denoise(:,iChan) = denoise;    ms_post_SW.deconv(:,iChan) = deconv;

    end
    toc;
    
    ms_pre_REM = MS_deconv2rate([], ms_pre_REM);
    ms_post_REM = MS_deconv2rate([], ms_post_REM);
    ms_pre_SW = MS_deconv2rate([], ms_pre_SW);
    ms_post_SW = MS_deconv2rate([], ms_post_SW);

            
    
    figure(1011)
    c_ord = linspecer(5);
    cla
    hold on
    bin_holder = zeros(length(all_binary_post_REM),1);
    tvec = (1:length(ms.denoise))/mode(diff(ms_seg_resize.time{1}));
    offset = .5;
    for ii = 1:20
        plot(tvec(logical(all_binary_post_REM(:,ii))), bin_holder(logical(all_binary_post_REM(:,ii))) +offset*ii, 'd','color',  c_ord(3,:), 'linewidth', 2)
        plot(tvec, all_RawTraces_post_REM(:, ii) +offset*ii, 'k', 'linewidth',1)
        plot(tvec, ms.denoise(:, ii) +offset*ii,'color',  c_ord(1,:), 'linewidth', 1)
        plot(tvec, ms.deconv(:, ii) +offset*ii,'color',  c_ord(2,:)', 'linewidth', 1)
        plot(tvec, ms.rate(:, ii) +offset*ii, 'color', c_ord(5,:))
    end
    
    %% convert Rate to Lui style firing intervals ( >3Sd)
    for ii = size(ms.rate,2):-1:1
        ms.rate_binary(:,ii) = ms.rate(:,ii) > (mean(ms.rate(:,ii)) + std(ms.rate(:,ii))*3);
    end
    
    %%
    
    load('all_binary_pre_SW.mat');
    load('all_binary_post_SW.mat');
    

    
    nFrames = 33; % number of frames before and after the center of each event.
    %% get the SCEs
    
    
    % get the NLX time for Ca data
    pre_tvec = []; post_tvec = [];
    for ii = 1:length(ms_seg_resize.NLX_evt)
        if strcmp(ms_seg_resize.pre_post{ii}, 'pre') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
            pre_tvec = [pre_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
        elseif strcmp(ms_seg_resize.pre_post{ii}, 'post') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
            post_tvec = [post_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
        end
    end
    % get the indicies for recodring blocks.
    pre_seg_idx = find(diff(pre_tvec) > 1)+1;
    post_seg_idx = find(diff(post_tvec) > 1)+1;
    
    frame_n = floor(nFrames/2);
    
    if mod(frame_n, 2)==0; frame_n = frame_n+1; end
    frame_n_200 = floor(nFrames/5);
    shuff = 100;
    
    % run for pre
    data_in = all_binary_pre_SW;
    tic
    all_shuff = [];
    for iS = shuff:-1:1
        this_data = [];
        for ii = size(data_in, 2):-1:1
            this_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
        end % end cells
        
        %     all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); %malvache
        all_shuff(iS, :) = ((movsum(nansum(this_data,1), frame_n_200))./size(data_in,2))*100; % liu
    end % end shuff
    toc
    
    % find time points in real data that exceed threshold
    %  liu pop % version.
    pop_act = ((movsum(nansum(data_in,2), frame_n_200))./size(data_in,2))*100; %
    shuff_95 = prctile(all_shuff, 95, 'all');
    
    thresh = shuff_95;
    % malvache
    % thresh = shuff_mean + 3*shuff_sd;
    % shuff_mean = mean(all_shuff,'all');
    % shuff_sd = std(all_shuff, [],'all');
    % pop_act = movmean(nansum(data_in,2), frame_n_200);
    
    
    % exlude events that are too close. use findpeaks
    [peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', nFrames);
    tvec_SCE_pre =  0:(1/nFrames):(length(data_in)-1)*(1/nFrames);
    if plot_f
        figure(101);
        clf
        ax1(1) =subplot(20,1,1:7);
        MS_Ca_Raster(data_in', tvec_SCE_pre);
        % xlabel('frame number')
        ylabel('cell id')
        set(gca, 'xtick', [])
        if ~isempty(pre_seg_idx)
            vline(tvec_SCE_pre(pre_seg_idx), '--w')
        end
        ax1(2) =subplot(20,1,8:9);
        hold on
        plot(tvec_SCE_pre, pop_act);
        plot(tvec_SCE_pre(SCE_idx), pop_act(SCE_idx), 'x')
        xlim([tvec_SCE_pre(1) tvec_SCE_pre(end)])
        hline(thresh)
        ylabel('% active cells')
        set(gca, 'xtick', [])
    end
    
    % collect pre
    pop_act_pre = pop_act;
    SCE_idx_pre = SCE_idx;
    shuff_pre = all_shuff;
    thresh_pre = thresh;
    
    
    % run for post
    data_in = all_binary_post_SW;
    tic
    all_shuff = [];
    for iS = shuff:-1:1
        this_data = [];
        for ii = size(data_in, 2):-1:1
            this_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
        end % end cells
        
        %    all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); %malvache
        all_shuff(iS, :) = ((movsum(nansum(this_data,1), frame_n_200))./size(data_in,2))*100; % liu
    end % end shuff
    toc
    
    % find time points in real data that exceed threshold
    %  liu pop % version.
    pop_act = ((movsum(nansum(data_in,2), frame_n_200))./size(data_in,2))*100; %
    shuff_95 = prctile(all_shuff, 95, 'all');
    
    thresh = shuff_95;
    % malvache
    % thresh = shuff_mean + 3*shuff_sd;
    % shuff_mean = mean(all_shuff,'all');
    % shuff_sd = std(all_shuff, [],'all');
    % pop_act = movmean(nansum(data_in,2), frame_n_200);
    
    % exlude events that are too close. use findpeaks
    [peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', nFrames);
    
    tvec_SCE_post =  0:(1/nFrames):(length(data_in)-1)*(1/nFrames);
    if plot_f
        ax2(1) =subplot(20,1,11:17);
        MS_Ca_Raster(data_in', tvec_SCE_post);
        % xlabel('frame number')
        ylabel('cell id')
        set(gca, 'xtick', [])
        if ~isempty(post_seg_idx)
            vline(tvec_SCE_post(post_seg_idx), '--w')
        end
        
        ax2(2) =subplot(20,1,18:19);
        hold on
        plot(tvec_SCE_post, pop_act);
        plot(tvec_SCE_post(SCE_idx), pop_act(SCE_idx), 'x')
        xlim([tvec_SCE_post(1) tvec_SCE_post(end)])
        hline(thresh)
        ylabel('% active cells')
        set(gca, 'xtick', [])
    end
    
    % collect post
    pop_act_post = pop_act;
    SCE_idx_post = SCE_idx;
    shuff_post = all_shuff;
    thresh_post = thresh;
    
    
    
    
    pre_SCE_times = pre_tvec(SCE_idx_pre);
    post_SCE_times = post_tvec(SCE_idx_post);
    
    
end