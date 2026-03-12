%% sandbox_CIHR_HF_2026_Recall

% load some intermediate data

load("HF1b3_TFC_D5.mat")


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


for ii = 1:length(data.vr.evt.t)
    data.vr.evt.t{ii} = (data.vr.evt.t{ii}/1000)+tfc_on.tstart(1); 
end

data.vr.pos.tvec = (data.vr.pos.tvec./1000)+tfc_on.tstart(1);
%% convert to evts  
labels = unique(data.vr.evt.label);

phases = {'baseline', 'tone1', 'tone2'};


log_iv.tone1 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')},...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone1')}+20); 

log_iv.tone2 = iv(data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')},...
    data.vr.evt.t{ismember(data.vr.evt.label, 'tone2')}+20); 

log_iv.baseline = iv(sort([log_iv.tone1.tstart; log_iv.tone2.tstart])-20,...
    sort([log_iv.tone1.tstart; log_iv.tone2.tstart])); 


log_iv.end = data.vr.evt.t{ismember(data.vr.evt.label, 'Finished')}; 


data.vr.pos.data(2,:) = diff([data.vr.pos.data data.vr.pos.data(end)]); %, {'encoder' 'movement'}); 
data.vr.pos.label = {'encoder' 'movement'}; 

% movement binary
mov_bin = data.vr.pos.data(2,:) > 0+realmin;
data.vr.pos.data(3,:) = mov_bin;
data.vr.pos.label{3} = 'move binary';

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
            trials.(['trace' phases{iP}(end)]).licks(ii) = length(this.t{1})./15;

            % movement
            this = restrict(vr_pos, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['trace' phases{iP}(end)]).mov(ii) = length(this.t{1})./15;

            % swr
            this = restrict(swr_times, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['trace' phases{iP}(end)]).swr(ii) = length(this.t{1})./15;


            % baselines
            this = restrict(vr_licks, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['baseline' phases{iP}(end)]).licks(ii) = length(this.t{1})./15;

            % movement
            this = restrict(vr_pos, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['baseline' phases{iP}(end)]).mov(ii) = length(this.t{1})./15;

            % swr
            this = restrict(swr_times, log_iv.(phases{iP}).tstart(ii)+20, log_iv.(phases{iP}).tend(ii)+15);
            trials.(['baseline' phases{iP}(end)]).swr(ii) = length(this.t{1})./15;

        end
    end
    
end


   

%% figure for the session to chekc everything
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

% tone 1
ax(1) = subplot(4,5,[1:3 6:8]);
cla
hold on
max_pos = max(data.vr.pos.data(2,:)); 

for ii = 1:length(log_iv.tone1.tstart)
    this_r = restrict(data.vr.pos, log_iv.tone1.tstart(ii)-20, log_iv.tone1.tend(ii)+20);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:)./max_pos, 'color', c_ord(ii,:));

    % add the licks as a raster
    this_l = restrict(vr_licks, log_iv.tone1.tstart(ii)-20, log_iv.tone1.tend(ii)+20);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end


xline([0 20 35])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(end-1,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

title("CS+ | context 2")

ylabel('Movement (a.u.)')

xlim([-20 56])
y_l = ylim;


% tone 2
ax(1) = subplot(4,5,[11:13 16:18]);
cla
hold on

for ii = 1:length(log_iv.tone2.tstart)
    this_r = restrict(data.vr.pos, log_iv.tone2.tstart(ii)-20, log_iv.tone2.tend(ii)+20);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:)./max_pos, 'color', c_ord(ii,:));

    % add the licks as a raster
    this_l = restrict(vr_licks, log_iv.tone2.tstart(ii)-20, log_iv.tone2.tend(ii)+20);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end


xline([0 20 35])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',s_ord(3,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',s_ord(end,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

ylabel('Movement (a.u.)')

xlim([-20 56])
y_l = ylim;
title("CS- | context 2")


 % bar plots summarizing movement and licks.
    subplot(4,5, [4]);
    hold on
    hb = MS_bar_w_err(trials.baseline1.mov, trials.tone1.mov, [.7 .7 .7; r_ord(2,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone1.mov, trials.trace1.mov, [r_ord(2,:); r_ord(end-1,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % ')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
    ylim([0 70])

    subplot(4,5,[14]);
    hold on
    hb = MS_bar_w_err(trials.baseline2.mov, trials.tone2.mov, [.7 .7 .7; s_ord(3,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone2.mov, trials.trace2.mov, [s_ord(3,:) ; s_ord(end,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % ')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
    ylim([0 70])

% swr
    subplot(4,5, [5]);
    hold on
    hb = MS_bar_w_err(trials.baseline1.swr, trials.tone1.swr, [.7 .7 .7; r_ord(2,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone1.swr, trials.trace1.swr, [r_ord(2,:); r_ord(end-1,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('SWR rate')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})

    subplot(4,5,[15 ]);
    hold on
    hb = MS_bar_w_err(trials.baseline2.swr, trials.tone2.swr, [.7 .7 .7; s_ord(3,:)], 1, 'ttest' ,[0 1]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(trials.tone2.swr, trials.trace2.swr, [s_ord(3,:) ; s_ord(end,:)], 1,'ttest' ,[1 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('SWR rate')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})




    %% compare Tone1 vs Tone 2

    figure(104)