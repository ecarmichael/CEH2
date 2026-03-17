%% sandbox_CIHR_HF_2026_Recall

% load some intermediate data

load("HF1b3_TFC_D5.mat")

data.vr = HF_load_VR('C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\TFC\H1b3\test\HF1b3_2026-03-05_18-34-55_D5_Test\arduino_log_1772754078.csv')
%% align the vr to the recording. 

tfc_on = iv(data.evts.t{ismember(data.evts.label, '5')}(1,:), data.evts.t{ismember(data.evts.label, '5')}(2,:)); 

% if the onset TTL from the arduino has more than one time, get the longest
if length(tfc_on.tstart) > 1
    ttl_on_len = tfc_on.tend - tfc_on.tstart; 
    % [~, keep_idx] = max(ttl_on_len); % get the longest block
    keep_idx = (max(ttl_on_len) == ttl_on_len);
    tfc_on.tstart(~keep_idx) = []; 
    tfc_on.tend(~keep_idx) = [];
end

if isnan(tfc_on.tend(1))
    tfc_on.tend(1) = data.csc.tvec(end); 
end
for ii = 1:length(data.vr.evt.t)
    data.vr.evt.t{ii} = (data.vr.evt.t{ii})+tfc_on.tstart(1); 
end

data.vr.pos.tvec = (data.vr.pos.tvec)+tfc_on.tstart(1);

data.vr.pos.data(2,:) = diff([data.vr.pos.data data.vr.pos.data(end)]); %, {'encoder' 'movement'}); 
data.vr.pos.label = {'encoder' 'movement'}; 

% movement binary
mov_bin = data.vr.pos.data(2,:) > 0+realmin;
data.vr.pos.data(3,:) = mov_bin;
data.vr.pos.label{3} = 'move binary';
%% convert the spikes to a rate
data.rate= MS_spike2rate(data.S, data.csc.tvec, .025, 0); 

data.rate = restrict(data.rate, tfc_on.tstart(1), tfc_on.tend(1));

% movement binary
mov_ts = ts({data.vr.pos.tvec(mov_bin)}, {'move_binary'});
mov_rate =  MS_spike2rate(mov_ts, data.rate.tvec, .025, 0);
move_idx =  mov_rate.data > 0; 

if length(move_idx) ~= length(data.rate.tvec)
    pad = length(move_idx) - length(data.rate.tvec);
    move_idx(end+(abs(pad))) = 0; 
end
% detect assemblies
    % code_dir = ['C:\Users\' usr '\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Git_repos\Dos-Santos Assembly ICA'];

data.asmbly = MS_asmbly_ephys(data.rate.data, move_idx);


%% convert to evts  
labels = unique(data.vr.evt.label);

phases = {'baseline', 'tone1', 'tone2', 'tone1_c1', 'tone1_c2', 'tone2_c1', 'tone2_c2'};

context_idx = [1 1 1 2 2 1 1 1 2 2 1 1 2 2 2]; 
tone_idx = [1,1,2,1,2,2,1,2,1,1,2,1,2,2,1]; 
ctx1_t1 = tone_idx(context_idx == 1 ) == 1; % get the tone1 + context 1 indicies
ctx1_t2 = tone_idx(context_idx == 1 ) == 2; % get the tone2 + context 1 indicies
ctx2_t1 = tone_idx(context_idx == 2 ) == 1; % get the tone1 + context 2 indicies
ctx2_t2 = tone_idx(context_idx == 2 ) == 2; % get the tone2 + context 2 indicies


log_iv.tone1 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')},...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')}+20); 

log_iv.tone2 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')},...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')}+20); 

log_iv.baseline = iv(sort([log_iv.tone1.tstart; log_iv.tone2.tstart])-20,...
    sort([log_iv.tone1.tstart; log_iv.tone2.tstart])); 


% context 2
log_iv.tone1_c1 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')}(ctx1_t1),...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')}(ctx1_t1)+20); 

log_iv.tone1_c2 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')}(ctx2_t1),...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')}(ctx2_t1)+20); 

log_iv.tone2_c1 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')}(ctx1_t2),...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')}(ctx1_t2)+20); 

log_iv.tone2_c2 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')}(ctx2_t2),...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')}(ctx2_t2)+20); 

log_iv.baseline = iv(sort([log_iv.tone1.tstart; log_iv.tone2.tstart])-20,...
    sort([log_iv.tone1.tstart; log_iv.tone2.tstart])); 



log_iv.end = data.vr.evt.t{ismember(data.vr.evt.label, 'Finished')}; 


% get the SWR times and convert to ts
swr_times = ts({IVcenters(data.swr.iv)}, {'swr center'});  % Extract SWR center times
%% count the events per trial
trials = [];
vr_pos = ts({data.vr.evt.t{ismember(data.vr.evt.label, 'pos')}}, {'move'}); 
vr_licks = ts({data.vr.evt.t{ismember(data.vr.evt.label, 'Lick')}}, {'licks'}); 

for iP = 1:length(phases)
    % baselines
    for ii = 1:length(log_iv.(phases{iP}).tstart)
        this = restrict(vr_licks, log_iv.(phases{iP}).tstart(ii), log_iv.(phases{iP}).tend(ii));
        trials.(phases{iP}).licks(ii) = length(this.t{1})./(mode(log_iv.(phases{iP}).tend- log_iv.(phases{iP}).tstart));

        % movement
        this = restrict(vr_pos, log_iv.(phases{iP}).tstart(ii), log_iv.(phases{iP}).tend(ii));
        trials.(phases{iP}).mov(ii) = length(this.t{1})./(mode(log_iv.(phases{iP}).tend- log_iv.(phases{iP}).tstart));

        % swr
        this = restrict(swr_times, log_iv.(phases{iP}).tstart(ii), log_iv.(phases{iP}).tend(ii));
        trials.(phases{iP}).swr(ii) = length(this.t{1})./(mode(log_iv.(phases{iP}).tend- log_iv.(phases{iP}).tstart));

        if contains(phases{iP}, 'tone')

            this = restrict(vr_licks, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['trace' phases{iP}(5:end)]).licks(ii) = length(this.t{1})./15;

            % movement
            this = restrict(vr_pos, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['trace' phases{iP}(5:end)]).mov(ii) = length(this.t{1})./15;

            % swr
            this = restrict(swr_times, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['trace' phases{iP}(5:end)]).swr(ii) = length(this.t{1})./15;

            % % MUA
            % this = restrict(data.mua, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            % trials.(['trace' phases{iP}(5:end)]).mua(ii) = length(this.t{1})./15;
            % 

            % baselines
            this = restrict(vr_licks, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['baseline' phases{iP}(5:end)]).licks(ii) = length(this.t{1})./15;

            % movement
            this = restrict(vr_pos, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['baseline' phases{iP}(5:end)]).mov(ii) = length(this.t{1})./15;

            % swr
            this = restrict(swr_times, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['baseline' phases{iP}(5:end)]).swr(ii) = length(this.t{1})./15;

        end
    end
    
end

%% get the 1d tuning curves for each spike for the tone 1 and trace

% for ii  =size(data.rate.data,1):-1:1
% 
% 
% 
% end

%% assembly figure
r_ord = winter(6);
c_ord = MS_linspecer(9);

 MS_asmbly_ephys_raster(data.S, data.rate.tvec,data.asmbly.A_temp, data.asmbly.A_proj, 1:size(data.asmbly.A_proj,1))

 ylim([0 50])
yline(data.asmbly.w_thresh, '--k')
plot(data.csc.tvec, (data.csc.data(1,:)/50)-4, 'k')
for ii = 1:length(log_iv.baseline.tstart)
    rectangle(gca, 'Position',[log_iv.baseline.tstart(ii) -1 log_iv.baseline.tstart(ii) - log_iv.baseline.tstart(ii) 1], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
end
for ii = 1:length(log_iv.tone1.tstart)
    rectangle(gca, 'Position',[log_iv.tone1.tstart(ii) -1 log_iv.tone1.tend(ii)-log_iv.tone1.tstart(ii) 1], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
        rectangle(gca, 'Position',[log_iv.tone1.tend(ii) -1 15 1], 'FaceColor',r_ord(3,:), 'EdgeColor','none' )

end
for ii = 1:length(log_iv.tone2.tstart)
    rectangle(gca, 'Position',[log_iv.tone2.tstart(ii) -1 log_iv.tone2.tend(ii)-log_iv.tone2.tstart(ii) 1], 'FaceColor',r_ord(4,:), 'EdgeColor','none' )
        rectangle(gca, 'Position',[log_iv.tone2.tend(ii) -1 15 1], 'FaceColor',r_ord(5,:), 'EdgeColor','none' )

end
%% figure for the session to check everything
r_ord = winter(6);
c_ord = MS_linspecer(9);

figure(999)
clf
hold on

plot(data.vr.pos.tvec, data.vr.pos.data(2,:)./max(data.vr.pos.data(2,:)), 'k')
plot(data.vr.pos.tvec, data.vr.pos.data(2,:)./max(data.vr.pos.data(2,:)), 'k')

plot([vr_licks.t{1} vr_licks.t{1}]', [zeros(size(vr_licks.t{1}))-.5 zeros(size(vr_licks.t{1}))]', 'Color','r', 'LineWidth',2)
plot(data.licks.tvec, (data.licks.data(1,:)./max(data.licks.data(1,:))./4)-.5, 'k')

plot(data.csc.tvec, (data.csc.data(1,:)./max(data.csc.data(1,:))./2)-.75, 'Color', c_ord(2,:))
plot([swr_times.t{1} swr_times.t{1}]', [zeros(size(swr_times.t{1}))-.75 zeros(size(swr_times.t{1}))-.5]', 'Color','c', 'LineWidth',2)

plot(data.mua.tvec, (data.mua.data(1,:)./max(data.mua.data(1,:))./2)-1.5, 'Color', c_ord(2,:))

for ii = 1:length(log_iv.baseline.tstart)
    rectangle(gca, 'Position',[log_iv.baseline.tstart(ii) -.25 log_iv.baseline.tstart(ii) - log_iv.baseline.tstart(ii) .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
end
for ii = 1:length(log_iv.tone1.tstart)
    rectangle(gca, 'Position',[log_iv.tone1.tstart(ii) -.25 log_iv.tone1.tend(ii)-log_iv.tone1.tstart(ii) .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
end
for ii = 1:length(log_iv.tone2.tstart)
    rectangle(gca, 'Position',[log_iv.tone2.tstart(ii) -.25 log_iv.tone2.tend(ii)-log_iv.tone2.tstart(ii) .25], 'FaceColor',r_ord(4,:), 'EdgeColor','none' )
end

xlim([tfc_on.tstart(1), log_iv.end])
set(gca, 'color', [1 1 1], 'xtick', [-.75 -.5 ])


%% trialified plots

c_ord = MS_linspecer(9);

r_ord = winter(6);
s_ord = inferno(6);

figure(103)

clf

% tone 1 context 1
ax(1) = subplot(4,5,[1:3]);
cla
hold on
max_pos = max(data.vr.pos.data(2,:)); 
max_mua = max(data.mua.data(1,:)); 

for ii = 1:length(log_iv.tone1_c1.tstart)
    this_r = restrict(data.vr.pos, log_iv.tone1_c1.tstart(ii)-20, log_iv.tone1_c1.tend(ii)+35);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:)./max_pos, 'color', c_ord(ii,:));

    this_r = restrict(data.mua, log_iv.tone1_c1.tstart(ii)-20, log_iv.tone1_c1.tend(ii)+35);
    plot(this_r.tvec - this_r.tvec(1)-20, (this_r.data(1,:)./max_mua)-ii-.5, 'color', c_ord(ii,:));

    % add the licks as a raster
    this_l = restrict(vr_licks, log_iv.tone1_c1.tstart(ii)-20, log_iv.tone1_c1.tend(ii)+35);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end


xline([0 20 35])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(end-1,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 15 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

title("CS+ | context 1")
ylabel('Movement (a.u.)')
xlim([-20 56])
y_l = ylim;


% tone 1 context 2
ax(1) = subplot(4,5,[6:8]);
cla
hold on
max_pos = max(data.vr.pos.data(2,:)); 

for ii = 1:length(log_iv.tone1_c2.tstart)
    this_r = restrict(data.vr.pos, log_iv.tone1_c2.tstart(ii)-20, log_iv.tone1_c2.tend(ii)+35);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:)./max_pos, 'color', c_ord(ii,:));

    % add the licks as a raster
    this_l = restrict(vr_licks, log_iv.tone1_c2.tstart(ii)-20, log_iv.tone1_c2.tend(ii)+35);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end


xline([0 20 35])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(end-1,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 15 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

title("CS+ | context 2")
ylabel('Movement (a.u.)')
xlim([-20 56])
y_l = ylim;

% tone 2
ax(1) = subplot(4,5,[11:13]);
cla
hold on

for ii = 1:length(log_iv.tone2_c1.tstart)
    this_r = restrict(data.vr.pos, log_iv.tone2_c1.tstart(ii)-20, log_iv.tone2_c1.tend(ii)+35);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:)./max_pos, 'color', c_ord(ii,:));

    % add the licks as a raster
    this_l = restrict(vr_licks, log_iv.tone2_c1.tstart(ii)-20, log_iv.tone2_c1.tend(ii)+35);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end


xline([0 20 35])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',s_ord(3,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',s_ord(end,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 15 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

ylabel('Movement (a.u.)')

xlim([-20 56])
y_l = ylim;
title("CS- | context 1")

% tone 2 context 2
ax(1) = subplot(4,5,[16:18]);
cla
hold on

for ii = 1:length(log_iv.tone2_c2.tstart)
    this_r = restrict(data.vr.pos, log_iv.tone2_c2.tstart(ii)-20, log_iv.tone2_c2.tend(ii)+35);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:)./max_pos, 'color', c_ord(ii,:));

    % add the licks as a raster
    this_l = restrict(vr_licks, log_iv.tone2_c2.tstart(ii)-20, log_iv.tone2_c2.tend(ii)+35);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end

xline([0 20 35])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',s_ord(3,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',s_ord(end,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 15 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

ylabel('Movement (a.u.)')

xlim([-20 56])
y_l = ylim;
title("CS- | context 2")


 % bar plots summarizing movement
    subplot(4,5, 4);
    hold on
    hb = MS_bar_w_err(trials.baseline1_c1.mov, trials.tone1_c1.mov, [.7 .7 .7; r_ord(2,:)], 1, 'ttest' ,[0 1]);

    hb = MS_bar_w_err(trials.tone1_c1.mov, trials.trace1_c1.mov, [r_ord(2,:); r_ord(end-1,:)], 1,'ttest' ,[1 2]);

    ylabel('mov % '); xlabel('Context 1')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline', 'CS+', 'CS+ trace'})
    ylim([0 70])

    % context 2
    subplot(4,5, 9);
    hold on
    hb = MS_bar_w_err(trials.baseline1_c2.mov, trials.tone1_c2.mov, [.7 .7 .7; r_ord(2,:)], 1, 'ttest' ,[0 1]);

    hb = MS_bar_w_err(trials.tone1_c2.mov, trials.trace1_c2.mov, [r_ord(2,:); r_ord(end-1,:)], 1,'ttest' ,[1 2]);

    ylabel('mov % '); xlabel('Context 2')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline', 'CS+', 'CS+ trace'})
    ylim([0 70])

% context 1
    subplot(4,5,14);
    hold on
    hb = MS_bar_w_err(trials.baseline2_c1.mov, trials.tone2_c1.mov, [.7 .7 .7; s_ord(3,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone2_c1.mov, trials.trace2_c1.mov, [s_ord(3,:) ; s_ord(end,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % '); xlabel('context 1')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
    ylim([0 70])

    % context 2
    subplot(4,5,19);
    hold on
    hb = MS_bar_w_err(trials.baseline2_c2.mov, trials.tone2_c2.mov, [.7 .7 .7; s_ord(3,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone2_c2.mov, trials.trace2_c2.mov, [s_ord(3,:) ; s_ord(end,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % '); xlabel('context 2')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
    ylim([0 70])

% swr
    subplot(4,5, 5);
    hold on
    hb = MS_bar_w_err(trials.baseline1_c1.swr, trials.tone1_c1.swr, [.7 .7 .7; r_ord(2,:)], 1, 'ttest' ,[0 1]);

    hb = MS_bar_w_err(trials.tone1_c1.swr, trials.trace1_c1.swr, [r_ord(2,:); r_ord(end-1,:)], 1,'ttest' ,[1 2]);

    ylabel('swr'); xlabel('Context 1')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline', 'CS+', 'CS+ trace'})
    ylim([0 .2])

    % context 2
    subplot(4,5, 10);
    hold on
    hb = MS_bar_w_err(trials.baseline1_c2.swr, trials.tone1_c2.swr, [.7 .7 .7; r_ord(2,:)], 1, 'ttest' ,[0 1]);

    hb = MS_bar_w_err(trials.tone1_c2.swr, trials.trace1_c2.swr, [r_ord(2,:); r_ord(end-1,:)], 1,'ttest' ,[1 2]);

    ylabel('swr'); xlabel('Context 2')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline', 'CS+', 'CS+ trace'})
    ylim([0 .2])

% context 1
    subplot(4,5,15);
    hold on
    hb = MS_bar_w_err(trials.baseline2_c1.swr, trials.tone2_c1.swr, [.7 .7 .7; s_ord(3,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone2_c1.swr, trials.trace2_c1.swr, [s_ord(3,:) ; s_ord(end,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('swr'); xlabel('context 1')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
    ylim([0 .2])

    % context 2
    subplot(4,5,20);
    hold on
    hb = MS_bar_w_err(trials.baseline2_c2.swr, trials.tone2_c2.swr, [.7 .7 .7; s_ord(3,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone2_c2.swr, trials.trace2_c2.swr, [s_ord(3,:) ; s_ord(end,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('swr'); xlabel('context 2')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
    ylim([0 .2])

    %% compare Tone1 vs Tone 2

    figure(104)