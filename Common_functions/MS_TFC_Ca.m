function TFC_out = MS_TFC_Ca(data_in, signal,hd,  plot_flag)


% starting TS for shock onset in each subject
meta.JKA_01.start_idx = 158; 
meta.JKA_01.start_ts = 10466/1000; 
meta.JKA_01.shock_idx = 4389; 
meta.JKA_01.shock_ts = 287887/1000; 


meta.JKA_04.shock_idx = 4308; 

meta.JKA_05.shock_idx = 4286; 
meta.JKA_05.shock_ts = 285631/1000; 

meta.JKA_07.shock_idx = 4276; 


%% init
if nargin < 1
    data_in = cd;
    signal = 'RawTraces';
    hd = []; 
    plot_flag = 0;
elseif nargin < 2
    hd = []; 
    signal = 'RawTraces';
    plot_flag = 0;
elseif nargin < 3
    hd = []; 
    plot_flag = 0;
    elseif nargin < 4
    plot_flag = 0;
end



TFC1.baseline = [0 240];
TFC1.tone = [240 260; 402 422; 564 584];
TFC1.trace = [260 280; 422 442; 584 604];
TFC1.shock = [280 282; 442 444; 604 606];
TFC1.ITI = [282 402; 444 564; 606 706];

TFC2.baseline = [0 240];
TFC2.tone = [240 300; 480 540; 720 780];
TFC2.trace = [300 320; 540 560; 780 800] ;
TFC2.ITI = [320 480; 560 720; 800 960];

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

parts = strsplit(ms.dirName, '/');
info.subject = parts{end-1};
info.session = parts{end-2};


if isfield(ms, 'keep_idx')
    if length(ms.keep_idx) == size(ms.denoise,2)
    ms = MS_Ca_good_cells(ms);
    fprintf('ms data has been curated using %.0f/%.0f good cells\n', sum(ms.keep_idx), length(ms.keep_idx))
    end
end


s_t = meta.([info.subject(1:3) '_' info.subject(end-1:end)]).shock_ts - 280; 

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
    ms_tfc = MS_deconv2rate(cfg_d, ms_tfc);
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

%%
tone_s_idx = nearest_idx(TFC1.tone(:,1),ms_tfc.time);
tone_e_idx = nearest_idx(TFC1.tone(:,2),ms_tfc.time);%tone_s_idx + (ceil((1/mode(diff(ms_tfc.time)))) * (TFC1.tone(1,2) - TFC1.tone(1,1)));

tone_diff = min(tone_e_idx - tone_s_idx);

tone_e_idx = tone_s_idx+ tone_diff;

trace_s_idx = nearest_idx(TFC1.trace(:,1),ms_tfc.time);
trace_e_idx = nearest_idx(TFC1.trace(:,2),ms_tfc.time);%tone_s_idx + (ceil((1/mode(diff(ms_tfc.time)))) * (TFC1.tone(1,2) - TFC1.tone(1,1)));

trace_diff = min(trace_e_idx - trace_s_idx);

trace_e_idx = trace_s_idx+ trace_diff;

shock_s_idx = nearest_idx(TFC1.shock(:,1),ms_tfc.time);
shock_e_idx = nearest_idx(TFC1.shock(:,2),ms_tfc.time);

shock_diff = min(shock_e_idx - shock_s_idx);

shock_e_idx = shock_s_idx+ shock_diff;

% tone mean avtivity
tone_mat = mean(cat(3,ms_tfc.(signal)(tone_s_idx(1):tone_e_idx(1),:),...
    ms_tfc.(signal)(tone_s_idx(2):tone_e_idx(2),:), ms_tfc.(signal)(tone_s_idx(3):tone_e_idx(3),:)),3)';

% tone_mat = tone_mat./max(tone_mat, [], 2)

% tracemean avtivity
trace_mat = mean(cat(3,ms_tfc.(signal)(trace_s_idx(1):trace_e_idx(1),:),...
    ms_tfc.(signal)(trace_s_idx(2):trace_e_idx(2),:), ms_tfc.(signal)(trace_s_idx(3):trace_e_idx(3),:)),3)';

shock_mat = mean(cat(3,ms_tfc.(signal)(shock_s_idx(1):shock_e_idx(1),:),...
    ms_tfc.(signal)(shock_s_idx(2):shock_e_idx(2),:), ms_tfc.(signal)(shock_s_idx(3):shock_e_idx(3),:)),3)';

% average motion for each period

if isfield(hd, 'motion')
    
    tone_motion = mean([hd.motion(tone_s_idx(1):tone_e_idx(1));hd.motion(tone_s_idx(2):tone_e_idx(2)); hd.motion(tone_s_idx(3):tone_e_idx(3))] ,1);
    trace_motion = mean([hd.motion(trace_s_idx(1):trace_e_idx(1));hd.motion(trace_s_idx(2):trace_e_idx(2)); hd.motion(trace_s_idx(3):trace_e_idx(3))] ,1);
    shock_motion = mean([hd.motion(shock_s_idx(1):shock_e_idx(1));hd.motion(shock_s_idx(2):shock_e_idx(2)); hd.motion(shock_s_idx(3):shock_e_idx(3))] ,1);
    
else

tone_motion = mean([hd.data(4,tone_s_idx(1):tone_e_idx(1));hd.data(4,tone_s_idx(2):tone_e_idx(2)); hd.data(4,tone_s_idx(3):tone_e_idx(3))] ,1); 
trace_motion = mean([hd.data(4,trace_s_idx(1):trace_e_idx(1));hd.data(4,trace_s_idx(2):trace_e_idx(2)); hd.data(4,trace_s_idx(3):trace_e_idx(3))] ,1); 
shock_motion = mean([hd.data(4,shock_s_idx(1):shock_e_idx(1));hd.data(4,shock_s_idx(2):shock_e_idx(2)); hd.data(4,shock_s_idx(3):shock_e_idx(3))] ,1); 
end
% sort the maps




[~, p] = max(trace_mat,[],2);

[~, s_idx]  = sort(p);

tone_mat_s = tone_mat(s_idx,:);

all_mat_s = ms_tfc.(signal)(:,s_idx);
trace_mat_s = trace_mat(s_idx,:);
shock_mat_s = shock_mat(s_idx,:);


% figure
% clf
% hold on
% imagesc(tone_mat_s)
% scatter(p(s_idx), 1:size(tone_mat,1), 'xr')


%%
if plot_flag
    clf
    ax(1) = subplot(5,3,1:3);
    cla
    hold on
    for ii = 1:size(TFC1.tone,1)
        rectangle('Position',[TFC1.tone(ii,1), min(ms_tfc.(signal)(:,1:50),[], 'all'),...
            TFC1.tone(ii,2) - TFC1.tone(ii,1) (max(ms_tfc.(signal)(:,1:50),[], 'all') - min(ms_tfc.(signal)(:,1:50),[], 'all'))],...
            'FaceColor',[c_ord(1,:) .6])
    end
    
    for ii = 1:size(TFC1.trace,1)
        rectangle('Position',[TFC1.trace(ii,1), min(ms_tfc.(signal)(:,1:50),[], 'all'),...
            TFC1.trace(ii,2) - TFC1.trace(ii,1) (max(ms_tfc.(signal)(:,1:50),[], 'all') - min(ms_tfc.(signal)(:,1:50),[], 'all'))],...
            'FaceColor',[c_ord(3,:) .6])
    end
    
    for ii = 1:size(TFC1.shock,1)
        rectangle('Position',[TFC1.shock(ii,1), min(ms_tfc.(signal)(:,1:50),[], 'all'),...
            TFC1.shock(ii,2) - TFC1.shock(ii,1) (max(ms_tfc.(signal)(:,1:50),[], 'all') - min(ms_tfc.(signal)(:,1:50),[], 'all'))],...
            'FaceColor',[c_ord(2,:) .6])
    end
    plot(ms_tfc.time, ms_tfc.(signal)(:,1:50), 'k')
        ylabel('dF/f')

    
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
    
    %%%%   plot the activity matrix
    ax(3) = subplot(5,3,7:9);
    cla
    hold on
    imagesc(ms_tfc.time,1:size(ms_tfc.(signal),2), all_mat_s' )
    ylabel('cell id')
            ylim([0 size(ms_tfc.(signal),2)]); 

    for ii = 1:size(TFC1.tone,1)
        rectangle('Position',[TFC1.tone(ii,1),0,...
            TFC1.tone(ii,2) - TFC1.tone(ii,1), size(ms_tfc.(signal),2)],...
            'EdgeColor',[c_ord(1,:) .6])
    end
    
    for ii = 1:size(TFC1.trace,1)
        rectangle('Position',[TFC1.trace(ii,1), 0,...
            TFC1.trace(ii,2) - TFC1.trace(ii,1) size(ms_tfc.(signal),2)],...
            'EdgeColor',[c_ord(3,:) .6])
    end
    
    for ii = 1:size(TFC1.shock,1)
        rectangle('Position',[TFC1.shock(ii,1), 0,...
            TFC1.shock(ii,2) - TFC1.shock(ii,1) size(ms_tfc.(signal),2)],...
            'EdgeColor',[c_ord(2,:) .6])
    end
%     xlim([0 ms_tfc.time(end)])
    %%%%% plot the tone trace and shock mats
    
    
    
    subplot(5,3,10)
    cla; hold on;
    imagesc(ms_tfc.time(tone_s_idx:tone_e_idx),  1:size(ms_tfc.(signal),2), tone_mat_s);
    xline(ms_tfc.time(tone_s_idx(1)))
    xline(ms_tfc.time(tone_e_idx(1)))
    xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])
        ylabel('cell id')
        ylim([0 size(ms_tfc.(signal),2)]); 

    subplot(5,3,11)
    cla; hold on;
    imagesc(ms_tfc.time(trace_s_idx:trace_e_idx), 1:size(ms_tfc.(signal),2), trace_mat_s);
    xline(ms_tfc.time(trace_s_idx(1)))
    xline(ms_tfc.time(trace_e_idx(1)))
    xlim([ms_tfc.time(trace_s_idx(1)) ms_tfc.time(trace_e_idx(1))])
            ylim([0 size(ms_tfc.(signal),2)]); 
        ylim([0 size(ms_tfc.(signal),2)]); 

    subplot(5,3,12)
    cla; hold on;
    imagesc(ms_tfc.time(shock_s_idx:shock_e_idx), 1:size(ms_tfc.(signal),2), shock_mat_s);
    xline(ms_tfc.time(shock_s_idx(1))+10)
    xline(ms_tfc.time(shock_e_idx(1))-10)
    xlim([ms_tfc.time(shock_s_idx(1)) ms_tfc.time(shock_e_idx(1))])
            ylim([0 size(ms_tfc.(signal),2)]); 

    subplot(5,3,13)
    cla; hold on;
    plot(ms_tfc.time(tone_s_idx:tone_e_idx), mean(tone_mat_s,1),'k', 'LineWidth',2);
    for ii = 1:3
        plot(ms_tfc.time(tone_s_idx:tone_e_idx), mean(ms_tfc.(signal)(tone_s_idx(ii):tone_e_idx(ii),:),2), 'LineWidth',1);
    end
    xline(ms_tfc.time(tone_s_idx(1))+10)
    xline(ms_tfc.time(tone_e_idx(1))-10)
    xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])
    ylabel('Mean activity')
    
    yyaxis right; cla
    plot(ms_tfc.time(tone_s_idx:tone_e_idx), tone_motion, 'LineWidth',2);
        legend({'mean', '1^{st}', '2^{nd}', '3^{rd}', 'motion'}, 'box', 'off', 'orientation', 'horizontal', 'location', 'north')

        
    subplot(5,3,14)
    cla; hold on;
    plot(ms_tfc.time(trace_s_idx:trace_e_idx),mean(trace_mat_s,1),'k', 'LineWidth',2);
    for ii = 1:3
        plot(ms_tfc.time(trace_s_idx:trace_e_idx), mean(ms_tfc.(signal)(trace_s_idx(ii):trace_e_idx(ii),:),2), 'LineWidth',1);
    end
    xline(ms_tfc.time(trace_s_idx(1))+10)
    xline(ms_tfc.time(trace_e_idx(1))-10)
    xlim([ms_tfc.time(trace_s_idx(1)) ms_tfc.time(trace_e_idx(1))])
    
    yyaxis right; cla
    plot(ms_tfc.time(trace_s_idx:trace_e_idx), trace_motion, 'LineWidth',2);
    
    subplot(5,3,15)
    cla; hold on;
    plot(ms_tfc.time(shock_s_idx:shock_e_idx), mean(shock_mat_s,1),'k', 'LineWidth',2);
    for ii = 1:3
        plot(ms_tfc.time(shock_s_idx:shock_e_idx), mean(ms_tfc.(signal)(shock_s_idx(ii):shock_e_idx(ii),:),2), 'LineWidth',1);
    end
    xlim([ms_tfc.time(shock_s_idx(1)) ms_tfc.time(shock_e_idx(1))])
    xline(ms_tfc.time(shock_s_idx(1))+10)
    xline(ms_tfc.time(shock_e_idx(1))-10)
    
    yyaxis right; cla
    plot(ms_tfc.time(shock_s_idx:shock_e_idx), shock_motion, 'LineWidth',2);
    ylabel('HD motion')
    
    
        ax(3) = subplot(5,3,7:9);
    linkaxes(ax, 'x')
        xlim([ms_tfc.time(1) ms_tfc.time(end)])

end

%% look for assemblies
ms_tfc.time = ms_tfc.time.*1000; 

%%
A = []; 

[A.A_temp, A.A_proj, A.data_h, A.tvec, A.opts, A.w_thresh, A.shuff_stats] = MS_PCA_ICA_only(ms_tfc, true(length(ms_tfc.time),1), 0.5);

[A.P_temp, A.P_proj, A.P_cells] = MS_Asmbly_select(A.A_temp, A.A_proj, 2)

%% plot 

MS_Asmbly_plot(A.P_temp, A.P_cells, [], [])

