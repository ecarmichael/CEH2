function Master_TFC_CA(data_dir)


%% init
if nargin < 1
    data_dir = cd; 
    signal = 'RawTraces'; 
    plot_flag = 0; 
end


cd(data_dir)

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
load('ms.mat')

% get the session info

parts = strsplit(ms.dirName, filesep); 
info.subject = parts{end-1}
info.session = parts{end-2}; 


if isfield(ms, 'keep_idx')
    ms = MS_Ca_good_cells(ms); 
fprintf('ms data has been curated using %.0f/%.0f good cells\n', sum(ms.keep_idx), length(ms.keep_idx))
end


%% split out the task and sleep


ms_pre =  MS_restrict(ms, ms.tvecs{1}(1),ms.tvecs{1}(end));

ms_tfc =  MS_restrict(ms, ms.tvecs{2}(1),ms.tvecs{2}(end));

ms_post =  MS_restrict(ms, ms.tvecs{3}(1),ms.tvecs{3}(end));

%%
if plot_flag 

    figure(101)

    subplot(3,3,1:3)
    hold on
    plot(ms_pre.time, ms_pre.(signal)(:,10), 'r')
        plot(ms_tfc.time, ms_tfc.(signal)(:,10), 'k')
            plot(ms_post.time, ms_post.(signal)(:,10), 'b')

end

%% re-zero the tfc

ms_tfc.time = ms_tfc.time- ms_tfc.time(1);

%%
tone_s_idx = nearest_idx(TFC1.tone(:,1)-10,ms_tfc.time); 
tone_e_idx = nearest_idx(TFC1.tone(:,2)+10,ms_tfc.time);%tone_s_idx + (ceil((1/mode(diff(ms_tfc.time)))) * (TFC1.tone(1,2) - TFC1.tone(1,1))); 

tone_diff = min(tone_e_idx - tone_s_idx); 

tone_e_idx = tone_s_idx+ tone_diff; 

trace_s_idx = nearest_idx(TFC1.trace(:,1)-10,ms_tfc.time); 
trace_e_idx = nearest_idx(TFC1.trace(:,2)+10,ms_tfc.time);%tone_s_idx + (ceil((1/mode(diff(ms_tfc.time)))) * (TFC1.tone(1,2) - TFC1.tone(1,1))); 

trace_diff = min(trace_e_idx - trace_s_idx); 

trace_e_idx = trace_s_idx+ trace_diff; 

shock_s_idx = nearest_idx(TFC1.shock(:,1)-10,ms_tfc.time); 
shock_e_idx = nearest_idx(TFC1.shock(:,2)+10,ms_tfc.time); 

shock_diff = min(shock_e_idx - shock_s_idx); 

shock_e_idx = shock_s_idx+ shock_diff; 

% tone mean avtivity
tone_mat = mean(cat(3,ms_tfc.(signal)(tone_s_idx(1):tone_e_idx(1),:),...
    ms_tfc.(signal)(tone_s_idx(2):tone_e_idx(2),:), ms_tfc.(signal)(tone_s_idx(3):tone_e_idx(3),:)),3); 

% tracemean avtivity
trace_mat = mean(cat(3,ms_tfc.(signal)(trace_s_idx(1):trace_e_idx(1),:),...
    ms_tfc.(signal)(trace_s_idx(2):trace_e_idx(2),:), ms_tfc.(signal)(trace_s_idx(3):trace_e_idx(3),:)),3); 

shock_mat = mean(cat(3,ms_tfc.(signal)(shock_s_idx(1):shock_e_idx(1),:),...
    ms_tfc.(signal)(shock_s_idx(2):shock_e_idx(2),:), ms_tfc.(signal)(shock_s_idx(3):shock_e_idx(3),:)),3); 


%%
if plot_flag

    subplot(4,3,4:6)
cla
hold on
    for ii = 1:size(TFC1.tone,1) 
    rectangle('Position',[TFC1.tone(ii,1), min(ms_tfc.(signal)(:,1:50),[], 'all'),...
        TFC1.tone(ii,2) - TFC1.tone(ii,1) (max(ms_tfc.(signal)(:,1:50),[], 'all') - min(ms_tfc.(signal)(:,1:50),[], 'all'))],...
        'FaceColor',c_ord(1,:), 'FaceAlpha',.6)
    end

    for ii = 1:size(TFC1.trace,1) 
    rectangle('Position',[TFC1.trace(ii,1), min(ms_tfc.(signal)(:,1:50),[], 'all'),...
        TFC1.trace(ii,2) - TFC1.trace(ii,1) (max(ms_tfc.(signal)(:,1:50),[], 'all') - min(ms_tfc.(signal)(:,1:50),[], 'all'))],...
        'FaceColor',c_ord(3,:), 'FaceAlpha',.6)
    end

        for ii = 1:size(TFC1.shock,1) 
    rectangle('Position',[TFC1.shock(ii,1), min(ms_tfc.(signal)(:,1:50),[], 'all'),...
        TFC1.shock(ii,2) - TFC1.shock(ii,1) (max(ms_tfc.(signal)(:,1:50),[], 'all') - min(ms_tfc.(signal)(:,1:50),[], 'all'))],...
        'FaceColor',c_ord(2,:), 'FaceAlpha',.6)
    end
        plot(ms_tfc.time, ms_tfc.(signal)(:,1:50), 'k')

        xlim([ms_tfc.time(1) ms_tfc.time(end)])
        subplot(4,3,7)
        cla; hold on;
        imagesc(ms_tfc.time(tone_s_idx:tone_e_idx),  1:size(ms_tfc.(signal),2), tone_mat);
        xline(ms_tfc.time(tone_s_idx(1))+10)
        xline(ms_tfc.time(tone_e_idx(1))-10)
        xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])

        subplot(4,3,8)
        cla; hold on;
        imagesc(ms_tfc.time(trace_s_idx:trace_e_idx), 1:size(ms_tfc.(signal),2), trace_mat);
        xline(ms_tfc.time(trace_s_idx(1))+10)
        xline(ms_tfc.time(trace_e_idx(1))-10)
        xlim([ms_tfc.time(trace_s_idx(1)) ms_tfc.time(trace_e_idx(1))])

        subplot(4,3,9)
        cla; hold on;
        imagesc(ms_tfc.time(shock_s_idx:shock_e_idx), 1:size(ms_tfc.(signal),2), shock_mat);
        xline(ms_tfc.time(shock_s_idx(1))+10)
        xline(ms_tfc.time(shock_e_idx(1))-10)
        xlim([ms_tfc.time(shock_s_idx(1)) ms_tfc.time(shock_e_idx(1))])

        subplot(4,3,10);
        cla; hold on;
        plot(ms_tfc.time(tone_s_idx:tone_e_idx), mean(tone_mat,2),'k', 'LineWidth',2);
        for ii = 1:3
            plot(ms_tfc.time(tone_s_idx:tone_e_idx), mean(ms_tfc.(signal)(tone_s_idx(ii):tone_e_idx(ii),:),2), 'LineWidth',1);
        end
        xline(ms_tfc.time(tone_s_idx(1))+10)
        xline(ms_tfc.time(tone_e_idx(1))-10)
        xlim([ms_tfc.time(tone_s_idx(1)) ms_tfc.time(tone_e_idx(1))])
        legend({'mean', '1^{st}', '2^nd', '3^rd'})

        subplot(4,3,11);
        cla; hold on;
        plot(ms_tfc.time(trace_s_idx:trace_e_idx),mean(trace_mat,2),'k', 'LineWidth',2);
        for ii = 1:3
            plot(ms_tfc.time(trace_s_idx:trace_e_idx), mean(ms_tfc.(signal)(trace_s_idx(ii):trace_e_idx(ii),:),2), 'LineWidth',1);
        end
        xline(ms_tfc.time(trace_s_idx(1))+10)
        xline(ms_tfc.time(trace_e_idx(1))-10)
        xlim([ms_tfc.time(trace_s_idx(1)) ms_tfc.time(trace_e_idx(1))])

        subplot(4,3,12)
        cla; hold on;
        plot(ms_tfc.time(shock_s_idx:shock_e_idx), mean(shock_mat,2),'k', 'LineWidth',2);
        for ii = 1:3
            plot(ms_tfc.time(shock_s_idx:shock_e_idx), mean(ms_tfc.(signal)(shock_s_idx(ii):shock_e_idx(ii),:),2), 'LineWidth',1);
        end
        xlim([ms_tfc.time(shock_s_idx(1)) ms_tfc.time(shock_e_idx(1))])
        xline(ms_tfc.time(shock_s_idx(1))+10)
        xline(ms_tfc.time(shock_e_idx(1))-10)


end