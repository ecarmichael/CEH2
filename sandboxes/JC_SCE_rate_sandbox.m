clear all
close all

ms_dir = {'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\7_15_2019_PV1060_LTD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\7_19_2019_PV1060_LTD5', ...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\11_26_2019_PV1060_HATSwitch', ...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_12_2019_PV1069_LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1192\4_17_2021_PV1192_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_14_2021_pv1252_LTD3',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_16_2021_pv1252_LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_24_2021_pv1252_HATDSwitch',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_13_2021_pv1254_LTD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_15_2021_pv1254_LTD3',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_17_2021_pv1254_LTD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_19_2021_pv1254_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_23_2021_pv1254_HATD5',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_25_2021_pv1254_HATDSwitch',...
    };


plot_f = 0; % toggle plots for speed.

cutoff_n_cell = 50;

method = 'decon';

for iSess = 13:length(ms_dir)
    %% load some data
    
    cd(ms_dir{iSess})
    if ~exist('ms_trk.mat')
        continue
    end
    parts = strsplit(cd, filesep);
    if contains(parts{end}, 'pv')
        parts{end} =   strrep(parts{end}, 'pv', 'PV');
    end
    sess = strsplit(parts{end}, parts{end-1});
    sess_id= [parts{end-1} '_' sess{2}(2:end)];
    
    load('ms_resize.mat');
    load('ms_trk.mat');
    
    if strcmpi(method, 'decon')
        fprintf('<strong>%s</strong>: using method ''<strong>%s</strong>''...\n', mfilename, method)
        
        cd('C:\Users\ecarm\Documents\GitHub\OASIS_matlab')
        oasis_setup;
        cd(ms_dir{iSess})
        
        
        load('all_RawTraces_pre_REM.mat');
        load('all_RawTraces_post_REM.mat');
        
        load('all_RawTraces_pre_SW.mat');
        load('all_RawTraces_post_SW.mat');
        
        % decon the traces
        fprintf('\n<strong>%s</strong>: deconvolving traces...\n', mfilename)
        tic;
        for iChan = size(all_RawTraces_post_REM,2):-1:1
            disp(num2str(iChan))
            [denoise,deconv] = deconvolveCa(all_RawTraces_pre_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            temp_pre_REM.denoise(:,iChan) = denoise;    temp_pre_REM.deconv(:,iChan) = deconv;
            
            [denoise,deconv] = deconvolveCa(all_RawTraces_post_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            temp_post_REM.denoise(:,iChan) = denoise;    temp_post_REM.deconv(:,iChan) = deconv;
            
            [denoise,deconv] = deconvolveCa(all_RawTraces_pre_SW(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            temp_pre_SW.denoise(:,iChan) = denoise;    temp_pre_SW.deconv(:,iChan) = deconv;
            
            [denoise,deconv] = deconvolveCa(all_RawTraces_post_SW(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            temp_post_SW.denoise(:,iChan) = denoise;    temp_post_SW.deconv(:,iChan) = deconv;
            
            [denoise,deconv] = deconvolveCa(ms_trk.RawTraces(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            temp_trk.denoise(:,iChan) = denoise;    temp_trk.deconv(:,iChan) = deconv;
        end
        toc;
        
        temp_pre_REM = MS_deconv2rate([], temp_pre_REM);
        temp_post_REM = MS_deconv2rate([], temp_post_REM);
        temp_pre_SW = MS_deconv2rate([], temp_pre_SW);
        temp_post_SW = MS_deconv2rate([], temp_post_SW);
        temp_trk = MS_deconv2rate([], temp_trk);
        
        
        %% convert Rate to Lui style firing intervals ( >3Sd)
        for ii = size(temp_pre_REM.rate,2):-1:1
            % pre REM
            temp_pre_REM.rate_binary(:,ii) = temp_pre_REM.rate(:,ii) > (mean(temp_pre_REM.rate(:,ii)) + std(temp_pre_REM.rate(:,ii))*3);
            % post REM
            temp_post_REM.rate_binary(:,ii) = temp_post_REM.rate(:,ii) > (mean(temp_post_REM.rate(:,ii)) + std(temp_post_REM.rate(:,ii))*3);
            % ms_pre_SW
            temp_pre_SW.rate_binary(:,ii) = temp_pre_SW.rate(:,ii) > (mean(temp_pre_SW.rate(:,ii)) + std(temp_pre_SW.rate(:,ii))*3);
            % pre REM
            temp_post_SW.rate_binary(:,ii) = temp_post_SW.rate(:,ii) > (mean(temp_post_SW.rate(:,ii)) + std(temp_post_SW.rate(:,ii))*3);
            % track
            temp_trk.rate_binary(:,ii) = temp_trk.rate(:,ii) > (mean(temp_trk.rate(:,ii)) + std(temp_trk.rate(:,ii))*3);
        end
        
        data_deconv = [];
        data_decon.cfg = {'deconvolveCa(ms.RaceTraces(:,iChan), ''foopsi'', ''ar2'', ''smin'', -2.5, ''optimize_pars'', true, ''optimize_b'', true)'};
        data_deconv.pre_REM = temp_pre_REM;
        data_deconv.post_REM = temp_post_REM;
        
        data_deconv.trk = temp_trk;
        
        data_deconv.pre_SW = temp_pre_SW;
        data_deconv.post_SW = temp_post_SW;
        
        save('Deconv_data.mat', 'data_deconv')
        
        % use common data form for
        data.pre_REM = temp_pre_REM.rate_binary;
        data.post_REM = temp_post_REM.rate_binary;
        
        data.pre_SW = temp_pre_SW.rate_binary;
        data.post_SW = temp_post_SW.rate_binary;
        
        data.trk = temp_trk.rate_binary;
        %%   plot for debugging.
        %     figure(1011)
        %     c_ord = linspecer(5);
        %     cla
        %     hold on
        %     bin_holder = zeros(length(all_binary_post_REM),1);
        %     tvec = (1:length(ms_post_REM.denoise))/round(1/mode(diff(ms.time)));
        %     offset = .5;
        %     for ii = 600:620
        %         plot(tvec(logical(all_binary_post_REM(:,ii))), bin_holder(logical(all_binary_post_REM(:,ii))) +offset*ii, 'd','color',  c_ord(3,:), 'linewidth', 2)
        %         plot(tvec, all_detrendRaw_post_REM(:, ii) +offset*ii, 'k', 'linewidth',1)
        %         plot(tvec, ms_post_REM.denoise(:, ii) +offset*ii,'color',  c_ord(1,:), 'linewidth', 1)
        %         plot(tvec, ms_post_REM.deconv(:, ii) +offset*ii,'color',  c_ord(2,:)', 'linewidth', 1)
        %         plot(tvec, ms_post_REM.rate(:, ii) +offset*ii, 'color', c_ord(5,:))
        %         plot(tvec(ms_post_REM.rate_binary(:,ii)), bin_holder(ms_post_REM.rate_binary(:,ii)) +offset*ii, 's','color', c_ord(4,:))
        %     end
        %% clean up
        clear temp_* data_deconv ms_p* all_*
        
    else
        load('all_binary_pre_REM.mat');
        load('all_binary_post_REM.mat');
        load('all_binary_pre_SW.mat');
        load('all_binary_post_SW.mat');
        
        data.pre_REM = all_binary_pre_REM;
        data.post_REM = all_binary_post_REM;
        
        data.pre_SW = all_binary_pre_SW;
        data.post_SW = all_binary_post_SW;
        
        data.trk = ms_trk.Binary;
        clear all_*
        
    end
    
    
    nFrames = round(1/mode(diff(ms_trk.time/1000))); % number of frames before and after the center of each event.
    %% get the SCEs for each data type
    types = fieldnames(data);
    
    thresh = 5; % use basic 5% cutoff.
    trk_start = datevec(ms_trk.time_labels);
    trk_end = ms_trk.time(end)/1000;
    for iT = 1:length(types)
        
        % get a time vector
        %
        if strcmpi(types{iT}, 'trk')
            tvec = ms_trk.time;
            e_time = tvec;
        elseif contains(types{iT}, 'pre') && contains(types{iT}, 'REM')
            tvec = []; e_time = [];
            for ii = 1:length(ms_seg_resize.NLX_evt)
                if strcmp(ms_seg_resize.pre_post{ii}, 'pre') && strcmp(ms_seg_resize.hypno_label{ii}, 'REM')
                    t_time = ms_seg_resize.NLX_evt{ii}.t{end};
                    t_time = t_time - t_time(1);
                    e_time = [e_time t_time + etime(datevec(ms_seg_resize.time_labels{ii}), trk_start)];
                    tvec = [tvec ms_seg_resize.NLX_evt{ii}.t{end}];
                end
            end
        elseif contains(types{iT}, 'pre') && contains(types{iT}, 'SW')
            tvec = []; e_time = [];
            for ii = 1:length(ms_seg_resize.NLX_evt)
                if strcmp(ms_seg_resize.pre_post{ii}, 'pre') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
                    t_time = ms_seg_resize.NLX_evt{ii}.t{end};
                    t_time = t_time - t_time(1);
                    e_time = [e_time t_time + etime(datevec(ms_seg_resize.time_labels{ii}), trk_start)];
                    tvec = [tvec ms_seg_resize.NLX_evt{ii}.t{end}];
                end
            end
            
        elseif contains(types{iT}, 'post') && contains(types{iT}, 'REM')
            tvec = []; e_time = [];
            for ii = 1:length(ms_seg_resize.NLX_evt)
                if strcmp(ms_seg_resize.pre_post{ii}, 'post') && strcmp(ms_seg_resize.hypno_label{ii}, 'REM')
                    t_time = ms_seg_resize.NLX_evt{ii}.t{end};
                    t_time = t_time - t_time(1);
                    e_time = [e_time t_time + etime(datevec(ms_seg_resize.time_labels{ii}), trk_start)+trk_end];
                    tvec = [tvec ms_seg_resize.NLX_evt{ii}.t{end}];
                end
            end
            
        elseif contains(types{iT}, 'post') && contains(types{iT}, 'SW')
            tvec = []; e_time = [];
            for ii = 1:length(ms_seg_resize.NLX_evt)
                if strcmp(ms_seg_resize.pre_post{ii}, 'post') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
                    t_time = ms_seg_resize.NLX_evt{ii}.t{end};
                    t_time = t_time - t_time(1);
                    e_time = [e_time t_time + etime(datevec(ms_seg_resize.time_labels{ii}), trk_start)+trk_end];
                    tvec = [tvec ms_seg_resize.NLX_evt{ii}.t{end}];
                end
            end
        end
        
        
        %     % get the NLX time for Ca data
        %     pre_tvec = []; post_tvec = [];
        %     for ii = 1:length(ms_seg_resize.NLX_evt)
        %         if strcmp(ms_seg_resize.pre_post{ii}, 'pre') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
        %             pre_tvec = [pre_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
        %         elseif strcmp(ms_seg_resize.pre_post{ii}, 'post') && strcmp(ms_seg_resize.hypno_label{ii}, 'SW')
        %             post_tvec = [post_tvec ms_seg_resize.NLX_evt{ii}.t{end}];
        %         end
        %     end
        %     % get the indicies for recodring blocks.
        %     pre_seg_idx = find(diff(pre_tvec) > 1)+1;
        %     post_seg_idx = find(diff(post_tvec) > 1)+1;
        
        frame_n = floor(nFrames/2);
        
        if mod(frame_n, 2)==0; frame_n = frame_n+1; end
        frame_n_200 = floor(nFrames/5);
        shuff = 100;
        
        % run for pre
        data_in = data.(types{iT});
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
        pop_act_nCell = (movsum(nansum(data_in,2), frame_n_200));
        shuff_95 = prctile(all_shuff, 95, 'all');
        
        thresh_95prct = shuff_95;
        
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
            
            ax1(2) =subplot(20,1,8:9);
            yyaxis right
            plot(tvec_SCE_pre, pop_act_nCell);
            xlim([tvec_SCE_pre(1) tvec_SCE_pre(end)])
            ylabel('N active cells')
            
            
            hold on
            yyaxis left
            plot(tvec_SCE_pre, pop_act);
            plot(tvec_SCE_pre(SCE_idx), pop_act(SCE_idx), 'x')
            ylabel('% active cells')
            
            hline(thresh)
            set(gca, 'xtick', [])
            linkaxes(ax1, 'x')
            
        end
        
        % collect pre
        pop_act_pre = pop_act;
        SCE_idx_pre = SCE_idx;
        %     shuff_pre = all_shuff;
        %     thresh_pre = thresh;
        
        pop_act_out.(sess_id).(types{iT}) = pop_act;
        SCE_times_out.(sess_id).(types{iT}) = e_time(SCE_idx);
        Shuff_out.(sess_id).(types{iT}) = all_shuff;
        SCE_rate_out.(sess_id).(types{iT}) = (length(SCE_idx) /length(tvec)/frame_n_200);
        
    end
    clearvars -except  ms_dir iSess  plot_f cutoff_n_cell method pop_act_out SCE_times_out Shuff_out SCE_rate_out
    
end
%% collect the oututs
types = fieldnames(SCE_rate_out.PV1060_HATSwitch);
sess_list = fieldnames(SCE_rate_out);
% all_SCE_rate_trk = nan(1, length(sess_list)); 
% all_SCE_rate_pre_SW = nan(1, length(sess_list)); 
% all_SCE_rate_pre_REM = nan(1, length(sess_list)); 
% all_SCE_rate_post_SW = nan(1, length(sess_list)); 
% all_SCE_rate_post_REM = nan(1, length(sess_list)); 
clear all_SCE*
for iS = length(sess_list):-1:1
    
    if isfield(SCE_rate_out.(sess_list{iS}), 'trk')
        all_SCE_rate_trk(iS) = SCE_rate_out.(sess_list{iS}).trk;
    else
        all_SCE_rate_trk(iS) = NaN;
    end
    
    if isfield(SCE_rate_out.(sess_list{iS}), 'pre_SW')
        all_SCE_rate_pre_SW(iS) = SCE_rate_out.(sess_list{iS}).pre_SW;
    else
        all_SCE_rate_pre_SW(iS) = NaN;
    end
    if isfield(SCE_rate_out.(sess_list{iS}), 'pre_REM')
        all_SCE_rate_pre_REM(iS) = SCE_rate_out.(sess_list{iS}).pre_REM;
    else
        all_SCE_rate_pre_REM(iS) = NaN;
    end
    if isfield(SCE_rate_out.(sess_list{iS}), 'post_SW')
        all_SCE_rate_post_SW(iS) = SCE_rate_out.(sess_list{iS}).post_SW;
    else
        all_SCE_rate_post_SW(iS) = NaN;
    end
    if isfield(SCE_rate_out.(sess_list{iS}), 'post_REM')
        all_SCE_rate_post_REM(iS) = SCE_rate_out.(sess_list{iS}).post_REM;
    else
        all_SCE_rate_post_REM(iS) = NaN;
    end
    
    
end



%% print a summary figure
b_in = ([all_SCE_rate_pre_REM; all_SCE_rate_pre_SW; all_SCE_rate_trk; all_SCE_rate_post_SW; all_SCE_rate_post_REM]'); 


figure(202)
hold on
bar(1:5, nanmean(b_in))
e_b = errorbar(1:5, nanmean(b_in),nanstd(b_in)./sqrt(length(b_in)));
e_b.LineStyle = 'None'; 


