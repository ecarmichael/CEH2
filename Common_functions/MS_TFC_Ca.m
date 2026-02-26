function TFC_out = MS_TFC_Ca(data_in, signal,hd, start_t, plot_flag)


% starting TS for shock onset in each subject. The time needs to come from
% the behaviour camera frame index [where the light turned off]. 
meta.JKA_01.start_idx = 158; 
meta.JKA_01.start_ts = 10466/1000; 
meta.JKA_01.shock_idx = 4389; 
meta.JKA_01.shock_ts = 287887/1000; 


meta.JKA_04.shock_idx = 4308; 

meta.JKA_05.shock_idx = 4286; 
meta.JKA_05.shock_ts = 285631/1000; 

meta.JKA_07.shock_idx = 4276; 

meta.CHNRA2_2000_TFC1.start_idx = 128; 
meta.CHNRA2_2000_TFC2.start_idx = 197; 



%% init
if nargin < 1
    data_in = cd;
    signal = 'RawTraces';
    hd = [];
    start_idx = [];
    plot_flag = 0;
elseif nargin < 2
    hd = [];
    signal = 'RawTraces';
    start_idx = [];
    plot_flag = 0;
elseif nargin < 3
    hd = [];
    plot_flag = 0;
    start_idx = [];
elseif nargin < 4
    start_idx = [];
    plot_flag = 0;
elseif nargin < 5
    plot_flag = 0;
end

% CHRNA_2000_shock_idx = [5264 8143 11023 13905 16784] %diff ~ 2880
% TFC1.baseline = [0 240];
% TFC1.tone = [240 260; 402 422; 564 584; 706 726; 848 868 ];
% TFC1.trace = [260 280; 422 442; 584 604; 726 746; 868 888];
% TFC1.shock = [280 282; 442 444; 604 606; 746 748; 888 890];
% TFC1.ITI = [282 402; 444 564; 606 706; 748 848; 890 980];

TFC1.baseline = [0 300]; 
TFC1.tone = [300 320; 492 512; 684 704; 876 896; 1068 1088]; 
TFC1.trace =    [320 340; 512 532; 704 724; 896 916; 1088 1108];
TFC1.shock =        [340 342; 532 534; 724 726; 916 918; 1108 1110];
TFC1.ITI =              [342 492; 534 684; 726 876; 918 1048; 1110  1240]; 

% TFC2.baseline = [0 240];
% TFC2.tone = [240 300; 480 540; 720 780];
% TFC2.trace = [300 320; 540 560; 780 800] ;
% TFC2.ITI = [320 480; 560 720; 800 960];

TFC2.baseline = [0 300]; 
TFC2.tone = [300 320; 490 510; 680 700; 870 890; 1060 1080]; 
TFC2.trace =    [320 340; 510 530; 700 720; 890 910; 1080 1100];
TFC2.ITI =              [340 490; 530 680; 720 870; 910 1040; 1100  1240]; 

TFC3.baseline = [0 300];

c_ord = MS_linspecer(4);





%% Load the data
if ~isstruct(data_in)
    cd(data_in)
    load('ms.mat')

    hd = MS_Load_v4HD(cd, 0);
elseif isfield(data_in, 'numNeurons')
    ms = data_in; 
elseif isfield(data_in, 'ms')
    ms = data_in.ms; 
    if isfield(data_in, 'hd')
        hd = data_in.hd;
    end
end
clear data_in; 

% get the session info

parts = strsplit(ms.dirName, filesep);
info.subject = parts{end-2};
info.session = parts{end-1};


if isfield(ms, 'keep_idx')
    if length(ms.keep_idx) == size(ms.denoise,2)
    ms = MS_Ca_good_cells(ms);
    fprintf('ms data has been curated using %.0f/%.0f good cells\n', sum(ms.keep_idx), length(ms.keep_idx))
    end
end

if isempty(start_t)
s_t = meta.([info.subject(1:3) '_' info.subject(end-1:end)]).shock_ts - 280; 
end
%% split out the task and sleep

if isfield(ms, 'tvecs') && size(ms.tvecs,2) > 1
    
    ms_pre =  MS_restrict(ms, ms.tvecs{1}(1),ms.tvecs{1}(end));
    
    ms_tfc =  MS_restrict(ms, ms.tvecs{2}(1),ms.tvecs{2}(end));
    
    
    ms_post =  MS_restrict(ms, ms.tvecs{3}(1),ms.tvecs{3}(end));
    
    if iscell(hd) 
        hd_pre =  hd{1}; 
        hd_post =  hd{3}; 

        hd =  hd{2}; 
        
        
    end
else
    
    ms_tfc = ms;
    cfg_d = [];
    cfg_d.bins = 1;
    cfg_d.gau_win = 2; % window size for gaussian conv
    cfg_d.gau_sd = 1; % standard devision for gaussian conv.
    % ms_tfc = MS_deconv2rate(cfg_d, ms_tfc);
end
% %
% if plot_flag
%
%     figure(101)
%
%     subplot(3,3,1:3)
%     hold on
%     plot(ms_pre.time, ms_pre.(signal)(:,10), 'r')
%         plot(ms_tfc.time, ms_tfc.(signal)(:,10), 'k')
%             plot(ms_post.time, ms_post.(signal)(:,10), 'b')
%
% end



%% re-zero the tfc

ms_tfc.time = ms_tfc.time- ms_tfc.time(1);
hd.tvec = hd.tvec - hd.tvec(1); 

hd.tvec = hd.tvec - start_t; 
ms_tfc.time  = ms_tfc.time-start_t; 

hd_r = restrict(hd, 0, hd.tvec(end)); 
ms_tfc_r = MS_restrict(ms_tfc, 0, ms_tfc.time(end));

%%
tone_s_idx = nearest_idx(TFC1.tone(:,1),ms_tfc_r.time);
tone_e_idx = nearest_idx(TFC1.tone(:,2),ms_tfc_r.time);%tone_s_idx + (ceil((1/mode(diff(ms_tfc.time)))) * (TFC1.tone(1,2) - TFC1.tone(1,1)));

tone_diff = min(tone_e_idx - tone_s_idx);

tone_e_idx = tone_s_idx+ tone_diff;

trace_s_idx = nearest_idx(TFC1.trace(:,1),ms_tfc_r.time);
trace_e_idx = nearest_idx(TFC1.trace(:,2),ms_tfc_r.time);%tone_s_idx + (ceil((1/mode(diff(ms_tfc.time)))) * (TFC1.tone(1,2) - TFC1.tone(1,1)));

trace_diff = min(trace_e_idx - trace_s_idx);

trace_e_idx = trace_s_idx+ trace_diff;

shock_s_idx = nearest_idx(TFC1.shock(:,1),ms_tfc_r.time);
shock_e_idx = nearest_idx(TFC1.shock(:,2),ms_tfc_r.time);

shock_diff = min(shock_e_idx - shock_s_idx);

shock_e_idx = shock_s_idx+ shock_diff;

% collect the activity across events
tone_mat = []; trace_mat = []; shock_mat = []; 
tone_motion = []; trace_motion = []; shock_motion = []; 

% define the signal to use
% this_signal = zscore(ms_tfc_r.(signal),[],1); 
this_signal = ms_tfc_r.(signal); 


for ii = 1:length(tone_s_idx)

    % tone mean avtivity
    tone_mat = cat(3,tone_mat, this_signal(tone_s_idx(ii):tone_e_idx(ii),:));

    % tracemean avtivity
    trace_mat = cat(3,trace_mat, this_signal(trace_s_idx(ii):trace_e_idx(ii),:));

    shock_mat = cat(3,shock_mat, this_signal(shock_s_idx(ii):shock_e_idx(ii),:));

    if isfield(hd_r, 'motion')
        tone_motion = [tone_motionl; hd_r.motion(tone_s_idx(ii):tone_e_idx(ii))];
        trace_motion = [trace_motion; hd_r.motion(trace_s_idx(ii):trace_e_idx(ii))];
        shock_motion = [shock_motion; hd_r.motion(shock_s_idx(ii):shock_e_idx(ii))];
    else
        tone_motion = [tone_motion; hd_r.data(4,tone_s_idx(1):tone_e_idx(1))];
        trace_motion = [trace_motion; hd_r.data(4,trace_s_idx(1):trace_e_idx(1))];
        shock_motion = [shock_motion; hd_r.data(4,shock_s_idx(1):shock_e_idx(1))];
    end
end

motion_max = max([tone_motion, trace_motion, shock_motion], [], 'all'); 

% tone_mat = tone_mat./max(tone_mat, [], 2)
% 
% % tracemean avtivity
% trace_mat = mean(cat(3,ms_tfc_r.(signal)(trace_s_idx(1):trace_e_idx(1),:),...
%     ms_tfc_r.(signal)(trace_s_idx(2):trace_e_idx(2),:), ms_tfc_r.(signal)(trace_s_idx(3):trace_e_idx(3),:),...
%     ms_tfc_r.(signal)(trace_s_idx(4):trace_e_idx(4),:),ms_tfc_r.(signal)(trace_s_idx(5):trace_e_idx(5),:)),3)';
% 
% shock_mat = mean(cat(3,ms_tfc_r.(signal)(shock_s_idx(1):shock_e_idx(1),:),...
%     ms_tfc_r.(signal)(shock_s_idx(2):shock_e_idx(2),:), ms_tfc_r.(signal)(shock_s_idx(3):shock_e_idx(3),:),...
%     ms_tfc_r.(signal)(shock_s_idx(4):shock_e_idx(4),:),ms_tfc_r.(signal)(shock_s_idx(5):shock_e_idx(5),:)),3)';
% 
% % average motion for each period


% % sort the maps
% [~, p] = max(trace_mat,[],2);
% 
% [~, s_idx]  = sort(p);
% 
% tone_mat_s = tone_mat(s_idx,:);
% 
% all_mat_s = ms_tfc_r.(signal);
% trace_mat_s = trace_mat(s_idx,:);
% shock_mat_s = shock_mat(s_idx,:);



%%
if plot_flag
    clf
    ax(1) = subplot(5,3,1:3);
    cla
    hold on
    for ii = 1:size(TFC1.tone,1)
        rectangle('Position',[TFC1.tone(ii,1), min(ms_tfc.(signal),[], 'all'),...
            TFC1.tone(ii,2) - TFC1.tone(ii,1) (max(ms_tfc.(signal),[], 'all') - min(ms_tfc.(signal),[], 'all'))],...
            'FaceColor',[c_ord(1,:) .6])
    end
    
    for ii = 1:size(TFC1.trace,1)
        rectangle('Position',[TFC1.trace(ii,1), min(ms_tfc.(signal),[], 'all'),...
            TFC1.trace(ii,2) - TFC1.trace(ii,1) (max(ms_tfc.(signal),[], 'all') - min(ms_tfc.(signal),[], 'all'))],...
            'FaceColor',[c_ord(3,:) .6])
    end
    
    for ii = 1:size(TFC1.shock,1)
        rectangle('Position',[TFC1.shock(ii,1), min(ms_tfc.(signal),[], 'all'),...
            TFC1.shock(ii,2) - TFC1.shock(ii,1) (max(ms_tfc.(signal),[], 'all') - min(ms_tfc.(signal),[], 'all'))],...
            'FaceColor',[c_ord(2,:) .6])
    end
    plot(ms_tfc.time, ms_tfc.(signal), 'k')
        ylabel('dF/f')
    xlim([0 ms_tfc.time(end)])
    set(gca, 'xtick', 0:300:max(ms_tfc.time), 'XTickLabel', 0:5:max(ms_tfc.time)/60)
    
    %$%% plot the head motion
    ax(2) = subplot(5,3,4:6);
    cla
    hold on
    plot(hd.tvec, hd.data(4,:), 'k')
    ylabel('Head motion')
    hd_m = max(hd.data(4,:))*1.1;  
    
    for ii = 1:size(TFC1.tone,1)
        rectangle('Position',[TFC1.tone(ii,1),0,...
            TFC1.tone(ii,2) - TFC1.tone(ii,1), hd_m],...
            'EdgeColor',[c_ord(1,:) .6])
    end
    
    for ii = 1:size(TFC1.trace,1)
        rectangle('Position',[TFC1.trace(ii,1), 0,...
            TFC1.trace(ii,2) - TFC1.trace(ii,1) hd_m],...
            'EdgeColor',[c_ord(3,:) .6])
    end
    
    for ii = 1:size(TFC1.shock,1)
        rectangle('Position',[TFC1.shock(ii,1), 0,...
            TFC1.shock(ii,2) - TFC1.shock(ii,1) hd_m],...
            'EdgeColor',[c_ord(2,:) .6])
    end

        xlim([0 hd.tvec(end)])

    %%%%   plot the activity matrix
    ax(3) = subplot(5,3,7:9);
    cla
    hold on
    imagesc(ms_tfc.time,1:size(ms_tfc.(signal),2), ms_tfc_r.(signal)' )
    ylabel('cell id')
            ylim([0.5 size(ms_tfc.(signal),2)]+.5); 

    for ii = 1:size(TFC1.tone,1)
        rectangle('Position',[TFC1.tone(ii,1),.5,...
            TFC1.tone(ii,2) - TFC1.tone(ii,1), size(ms_tfc.(signal),2)+.5],...
            'EdgeColor',[c_ord(1,:) .6])
    end
    
    for ii = 1:size(TFC1.trace,1)
        rectangle('Position',[TFC1.trace(ii,1), .5,...
            TFC1.trace(ii,2) - TFC1.trace(ii,1) size(ms_tfc.(signal),2)+.5],...
            'EdgeColor',[c_ord(3,:) .6])
    end
    
    for ii = 1:size(TFC1.shock,1)
        rectangle('Position',[TFC1.shock(ii,1), .5,...
            TFC1.shock(ii,2) - TFC1.shock(ii,1) size(ms_tfc.(signal),2)+.5],...
            'EdgeColor',[c_ord(2,:) .6])
    end
    xlim([0 ms_tfc.time(end)])
    
    %%%%% plot the tone trace and shock mats
   
    
    subplot(5,3,10)
    cla; hold on;
    this_t = ms_tfc.time(tone_s_idx(1):tone_e_idx(1)); 
    this_t = this_t-this_t(1); 

    imagesc(this_t,  1:size(ms_tfc.(signal),2), mean(tone_mat,3)');
    xlim([this_t(1) ceil(this_t(end))])
        ylabel('cell id')
        ylim([0+.5 size(ms_tfc.(signal),2)+.5]); 
        yyaxis right
        plot(this_t, mean(tone_motion), '-w')
        ylabel('mean motion')
        ylim([0 motion_max])


    subplot(5,3,13)
    cla; hold on;
    for ii = 1:length(tone_s_idx)
        plot(this_t, mean(ms_tfc.(signal)(tone_s_idx(ii):tone_e_idx(ii),:),2), 'LineWidth',.5);
    end
        plot(this_t, mean(mean(tone_mat,3), 2)','k', 'LineWidth',3);

    % xline(ms_tfc.time(tone_s_idx(1))+10)
    % xline(ms_tfc.time(tone_e_idx(1))-10)
    % xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])
    ylabel('Mean activity')
    
    yyaxis right; cla
    plot(this_t, tone_motion,'-r', 'LineWidth',2);
        legend({'1^{st}', '2^{nd}', '3^{rd}', '4^{rd}', '5^{rd}','mean', 'motion'}, 'box', 'off', 'orientation', 'horizontal', 'location', 'north')
    set(gca, 'ycolor', 'r')
title('tone')


    subplot(5,3,11)
    cla; hold on;
    this_t = ms_tfc.time(trace_s_idx(1):trace_e_idx(1)); 
    this_t = this_t-this_t(1); 

    imagesc(this_t,  1:size(ms_tfc.(signal),2), mean(trace_mat,3)');
    xlim([this_t(1) ceil(this_t(end))])
        ylabel('cell id')
        ylim([0+.5 size(ms_tfc.(signal),2)+.5]); 
      yyaxis right
        plot(this_t,mean(trace_motion), '-w')
        ylabel('mean motion')
        ylim([0 motion_max])


    subplot(5,3,14)
    cla; hold on;
    for ii = 1:length(tone_s_idx)
        plot(this_t, mean(ms_tfc.(signal)(trace_s_idx(ii):trace_e_idx(ii),:),2), 'LineWidth',.5);
    end
        plot(this_t, mean(mean(trace_mat,3), 2)','k', 'LineWidth',3);

    % xline(ms_tfc.time(tone_s_idx(1))+10)
    % xline(ms_tfc.time(tone_e_idx(1))-10)
    % xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])
    ylabel('Mean activity')
    
    yyaxis right; cla
    plot(this_t, trace_motion,'-r', 'LineWidth',2);
        legend({'1^{st}', '2^{nd}', '3^{rd}', '4^{rd}', '5^{rd}','mean', 'motion'}, 'box', 'off', 'orientation', 'horizontal', 'location', 'north')
    set(gca, 'ycolor', 'r')
title('trace')

    subplot(5,3,12)
    cla; hold on;
    this_t = ms_tfc.time(shock_s_idx(1):shock_e_idx(1)); 
    this_t = this_t-this_t(1); 

    imagesc(this_t,  1:size(ms_tfc.(signal),2), mean(shock_mat,3)');
    xlim([this_t(1) ceil(this_t(end))])
        ylabel('cell id')
        ylim([0+.5 size(ms_tfc.(signal),2)+.5]); 
      yyaxis right
        plot(this_t,mean(shock_motion), '-w')
        ylabel('mean motion')
        ylim([0 motion_max])

    subplot(5,3,15)
    cla; hold on;
    for ii = 1:length(tone_s_idx)
        plot(this_t, mean(ms_tfc.(signal)(shock_s_idx(ii):shock_e_idx(ii),:),2), 'LineWidth',.5);
    end
        plot(this_t, mean(mean(shock_mat,3), 2)','k', 'LineWidth',3);

    ylabel('Mean activity')
    yyaxis right; cla
    plot(this_t, shock_motion,'-r', 'LineWidth',2);
        legend({'1^{st}', '2^{nd}', '3^{rd}', '4^{rd}', '5^{rd}','mean', 'motion'}, 'box', 'off', 'orientation', 'horizontal', 'location', 'north')
    set(gca, 'ycolor', 'r')
    title('shock')
    % 
    % 
    % subplot(5,3,13)
    % cla; hold on;
    % plot(mean(ms_tfc.time(tone_s_idx:tone_e_idx)), mean(tone_mat_s,1),'k', 'LineWidth',2);
    % for ii = 1:3
    %     plot(ms_tfc.time(tone_s_idx:tone_e_idx), mean(ms_tfc.(signal)(tone_s_idx(ii):tone_e_idx(ii),:),2), 'LineWidth',1);
    % end
    % xline(ms_tfc.time(tone_s_idx(1))+10)
    % xline(ms_tfc.time(tone_e_idx(1))-10)
    % xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])
    % ylabel('Mean activity')
    % 
    % yyaxis right; cla
    % plot(ms_tfc.time(tone_s_idx:tone_e_idx), tone_motion, 'LineWidth',2);
    %     legend({'mean', '1^{st}', '2^{nd}', '3^{rd}', 'motion'}, 'box', 'off', 'orientation', 'horizontal', 'location', 'north')
    % 
    % 
    % subplot(5,3,14)
    % cla; hold on;
    % plot(ms_tfc.time(trace_s_idx:trace_e_idx),mean(trace_mat_s,1),'k', 'LineWidth',2);
    % for ii = 1:3
    %     plot(ms_tfc.time(trace_s_idx:trace_e_idx), mean(ms_tfc.(signal)(trace_s_idx(ii):trace_e_idx(ii),:),2), 'LineWidth',1);
    % end
    % xline(ms_tfc.time(trace_s_idx(1))+10)
    % xline(ms_tfc.time(trace_e_idx(1))-10)
    % xlim([ms_tfc.time(trace_s_idx(1)) ms_tfc.time(trace_e_idx(1))])
    % 
    % yyaxis right; cla
    % plot(ms_tfc.time(trace_s_idx:trace_e_idx), trace_motion, 'LineWidth',2);
    % 
    % subplot(5,3,15)
    % cla; hold on;
    % plot(ms_tfc.time(shock_s_idx:shock_e_idx), mean(shock_mat_s,1),'k', 'LineWidth',2);
    % for ii = 1:3
    %     plot(ms_tfc.time(shock_s_idx:shock_e_idx), mean(ms_tfc.(signal)(shock_s_idx(ii):shock_e_idx(ii),:),2), 'LineWidth',1);
    % end
    % xlim([ms_tfc.time(shock_s_idx(1)) ms_tfc.time(shock_e_idx(1))])
    % xline(ms_tfc.time(shock_s_idx(1))+10)
    % xline(ms_tfc.time(shock_e_idx(1))-10)
    % 
    % yyaxis right; cla
    % plot(ms_tfc.time(shock_s_idx:shock_e_idx), shock_motion, 'LineWidth',2);
    % ylabel('HD motion')
    
    
        ax(3) = subplot(5,3,7:9);
         linkaxes(ax, 'x')
        xlim([ms_tfc.time(1) ms_tfc.time(end)])

end

%% look for assemblies
% ms_tfc.time = ms_tfc.time.*1000; 

%%
% A = []; 
% 
% [A.A_temp, A.A_proj, A.data_h, A.tvec, A.opts, A.w_thresh, A.shuff_stats] = MS_PCA_ICA_only(ms_tfc, true(length(ms_tfc.time),1), 0.5);
% 
% [A.P_temp, A.P_proj, A.P_cells] = MS_Asmbly_select(A.A_temp, A.A_proj, 2)

%% plot 

% MS_Asmbly_plot(A.P_temp, A.P_cells, [], [])

