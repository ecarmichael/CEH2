%% sandbox_CIHR_HR_TFC

if ispc
    kilo_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-03_16-34-57_TFC_ENC\Record Node 113\experiment1\recording1\continuous\Intan_RHD_USB-100.Rhythm Data';
    OE_evts_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-03_16-34-57_TFC_ENC\Record Node 113\experiment1\recording1\events\Intan_RHD_USB-100.Rhythm Data\TTL';
    ctrl_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\conditioning_logs\';
    eye_dir = [];
elseif ismac
    kilo_dir = [];
    ctrl_dir ='/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/conditioning_logs/';
    ctrl_dir =  '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/conditioning_logs/';
    eye_dir = [];
end


%% load the control log

log_tab = readtable([ctrl_dir 'arduino_log_1756923023.csv']);

log_tab.time = log_tab.time - log_tab.time(1);  % zero to the first time point.

log_tab.time = log_tab.time/ 1000; % convert from ms to seconds;

% convert to evts

labels = unique(log_tab.phase);
phases = {'baseline', 'Tone', 'Trace', 'Puff'};

for ii = 1:length(phases)
    log.(phases{ii}) = log_tab.time(contains(log_tab.event, '-') & contains(log_tab.phase, phases{ii}));
end
log.end = log_tab.time(contains(log_tab.event, 'Finish'));

all_phases.t = log_tab.time(contains(log_tab.event, '-'));
all_phases.data = log_tab.phase(contains(log_tab.event, '-'));

% make intervals for each phase state
log_iv =[];
for ii = 1:length(phases)
    idx = find(ismember(all_phases.t, log.(phases{ii})));
    if max(idx) == length(all_phases.t)
        log_iv.(phases{ii}) = iv(all_phases.t(idx), [all_phases.t(idx(1:end-1))+1; log.end]);
    else
        log_iv.(phases{ii}) = iv(all_phases.t(idx), all_phases.t(idx+1));
    end
end

% get the encoder changes as a vector.

wheel.tvec = log_tab.time(contains(log_tab.event, 'pos'));
wheel.data = log_tab.encoderCount(contains(log_tab.event, 'pos'));

% interp
tvec_i = 0:.01:wheel.tvec(end);

wheel.data = interp1(wheel.tvec, wheel.data, tvec_i);
wheel.data(isnan(wheel.data)) = 0;
wheel.tvec = tvec_i;


wheel_tsd = tsd(wheel.tvec, [wheel.data; diff([wheel.data wheel.data(end)])], {'encoder' 'movement'}); 


% licks
licks = []; licks.type = 'ts';
licks.t{1} = log_tab.time(contains(log_tab.event, 'Lick'));
licks.label{1} = 'licks';
licks.cfg.history.mfun = {'ts'};
licks.cfg.history.cfg = {[]};

% movement binary

mov_bin = wheel.data(1,:) > 0.1;
wheel_tsd.data(3,:) = mov_bin;
wheel_tsd.label{3} = 'move binary';

% count the licks per trial
trials = [];
for ii = 1:length(log_iv.Tone.tstart)

    bl = restrict(licks, log_iv.Tone.tstart(ii)-20, log_iv.Tone.tstart(ii));
    trials.base.licks(ii) = length(bl.t{1})./20;

    tl = restrict(licks, log_iv.Tone.tstart(ii), log_iv.Tone.tend(ii));
    trials.tone.licks(ii) = length(tl.t{1})./20;

    trl = restrict(licks, log_iv.Trace.tstart(ii), log_iv.Trace.tend(ii));
    trials.trace.licks(ii) = length(trl.t{1})./15;

    % same thing but movement
    bl = restrict(wheel_tsd, log_iv.Tone.tstart(ii)-20, log_iv.Tone.tstart(ii));
    trials.base.mov(ii) = sum(bl.data(3,:)) / length(bl.tvec);

    tl = restrict(wheel_tsd, log_iv.Tone.tstart(ii), log_iv.Tone.tend(ii));
    trials.tone.mov(ii) = sum(tl.data(3,:)) / length(tl.tvec);

    trl = restrict(wheel_tsd, log_iv.Trace.tstart(ii), log_iv.Trace.tend(ii));
    trials.trace.mov(ii) = sum(trl.data(3,:)) / length(trl.tvec);

end

%% ephys
cd(kilo_dir)
OE_rec = readNPY('timestamps.npy');
OE_rec = [OE_rec(1) OE_rec(end)];

evts = OE_load_binary_evts(OE_evts_dir);

rec_iv = iv(evts.t{1}(1,:)-OE_rec(1), evts.t{1}(2,:)-OE_rec(1));
tone_iv = iv(evts.t{2}(1,:)-OE_rec(1), evts.t{2}(2,:)-OE_rec(1));
trace_iv = iv(evts.t{3}(1,:)-OE_rec(1), evts.t{3}(2,:)-OE_rec(1));


% get the good cells

S = OE_phy2TS(kilo_dir);

% restrict to only when the behaviour started.
S_r = restrict(S, rec_iv);

for ii = 1:length(S_r.t)
    S_r.t{ii} = S_r.t{ii} - rec_iv.tstart;
end

% figure
% plot(S_r)
% xline()

%% plot??
c_ord = MS_linspecer(8);

figure(101)

clf

ax(1) = subplot(5,1,1);
cla

plot(wheel_tsd.tvec, wheel_tsd.data(1,:), 'k');

for ii = 1:length(phases)
    xline(log.(phases{ii}), 'Color',c_ord(ii,:), 'LineWidth',3)
end
ylabel('Encoder pos')

ax(2) = subplot(5,1,2);
cla
plot(wheel_tsd.tvec, abs(wheel_tsd.data(2,:)), 'r');
ylabel('Movement (a.u.)')

for ii = 1:length(phases)
    xline(log.(phases{ii}), 'Color',c_ord(ii,:), 'LineWidth',3)
end


ax(3) = subplot(5,1,3);
cla
plot([licks.t{1} licks.t{1}]', [zeros(size(licks.t{1})) ones(size(licks.t{1}))]', 'k');
ylabel('Licks'); ylim([0 2]);

for ii = 1:length(phases)
    xline(log.(phases{ii}), 'Color',c_ord(ii,:), 'LineWidth',3)
end

ax(4) = subplot(5,1,4:5);
cla
plot(S_r);
% ylabel('Licks'); ylim([0 2]);

for ii = 1:length(phases)
    xline(log.(phases{ii}), 'Color',c_ord(ii,:), 'LineWidth',3)
end

linkaxes(ax, 'x')
xlim([0 log.end])


%% trial-ified the data

% wheel and movement

mov_bin = wheel.data(1,:) > 0.1; 
wheel.data(3,:) = mov_bin; 
wheel.label{3} = 'move binary';




%%
% plot??
c_ord = MS_linspecer(9);

r_ord = winter(6);
figure(102)

clf

ax(1) = subplot(3,2,[1 3]);
cla
hold on
for ii = 1:length(log_iv.Tone.tstart)
    this_r = restrict(wheel_tsd, log_iv.Tone.tstart(ii)-20, log_iv.Puff.tend(ii)+20);

    plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:), 'color', c_ord(ii,:));

    % add the licks as a raster

    this_l = restrict(licks, log_iv.Tone.tstart(ii)-20, log_iv.Puff.tend(ii)+20);
    plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
end


xline([0 20 35 36])
rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(4,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[35 -.25 1 .25], 'FaceColor',r_ord(6,:), 'EdgeColor','none' )
rectangle(gca, 'Position',[36 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

ylabel('Movement (a.u.)')

xlim([-20 56])
y_l = ylim;
% ylim([-.55 y_l(2)]);


% bar plots summarizing movement and licks.
ax(2) = subplot(3,2, 2);
hold on

hb = MS_bar_w_err(trials.tone.mov, trials.base.mov, c_ord(2:3,:), 1, 'ttest',[1 0]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

hb = MS_bar_w_err(trials.tone.mov, trials.trace.mov, c_ord(1:2,:), 1, 'ttest',[1 2]);

hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

ylabel('mov % ')

set(gca, 'XTick', 1:3, 'XTickLabel', {'Baseline [-20:0]', 'Tone', 'Trace'})

% licks mean.
ax(2) = subplot(3,2, 4);
hold on

hb = MS_bar_w_err(trials.tone.licks, trials.base.licks, c_ord(2:3,:), 1, 'ttest',[1 0]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

hb = MS_bar_w_err(trials.tone.licks, trials.trace.licks, c_ord(1:2,:), 1, 'ttest',[1 2]);

hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

linkaxes(ax, 'x')
xlim([0 log.end])
ylabel('mean licks rate')

set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'Tone', 'Trace'})


%% over Conditioning sessions

s_list = dir([ctrl_dir filesep 'arduino_log*Con*.csv']);


% loop
data_out = [];
for iS  = length(s_list):-1:1
    log_tab = readtable([ctrl_dir s_list(iS).name]);


    log_tab.time = log_tab.time - log_tab.time(1);  % zero to the first time point.
    log_tab.time = log_tab.time/ 1000; % convert from ms to seconds;

    % convert to evts

    labels = unique(log_tab.phase);
    if contains(s_list(iS).name, 'Con')
        phases = {'baseline', 'Tone', 'Trace', 'Puff'};
    else
        phases = {'baseline', 'tone1', 'tone2'};
    end

    for ii = 1:length(phases)
        log.(phases{ii}) = log_tab.time(contains(log_tab.event, '-') & contains(log_tab.phase, phases{ii}));
    end

    log.end = log_tab.time(contains(log_tab.event, 'Finish'));

    all_phases.t = log_tab.time(contains(log_tab.event, '-'));
    all_phases.data = log_tab.phase(contains(log_tab.event, '-'));

    % make intervals for each phase state
    log_iv =[];
    for ii = 1:length(phases)
        idx = find(ismember(all_phases.t, log.(phases{ii})));
        if max(idx) == length(all_phases.t)
            log_iv.(phases{ii}) = iv(all_phases.t(idx), [all_phases.t(idx(1:end-1))+1; log.end]);
        else
            log_iv.(phases{ii}) = iv(all_phases.t(idx), all_phases.t(idx+1));
        end
    end

    % get the encoder changes as a vector.

    wheel.tvec = log_tab.time(contains(log_tab.event, 'pos'));
    wheel.data = log_tab.encoderCount(contains(log_tab.event, 'pos'));

    % interp
    tvec_i = 0:.01:log.end;
    wheel.data = interp1(wheel.tvec, wheel.data, tvec_i);
    wheel.data(isnan(wheel.data)) = 0;
    wheel.tvec = tvec_i;

    % remove extreme values
    wheel_tsd.data(abs(wheel.data) > 10 ) = NaN;

    wheel.data = fillmissing(wheel.data, 'nearest');

    wheel_tsd = tsd(wheel.tvec, [wheel.data; diff([wheel.data wheel.data(end)])], {'encoder' 'movement'});

    % licks
    licks = []; licks.type = 'ts';
    licks.t{1} = log_tab.time(contains(log_tab.event, 'Lick'));
    licks.label{1} = 'licks';
    licks.cfg.history.mfun = {'ts'};
    licks.cfg.history.cfg = {[]};


    mov_bin = wheel.data(1,:) > 0.1;
    wheel_tsd.data(3,:) = mov_bin;
    wheel_tsd.label{3} = 'move binary';

    % count the licks per trial
    trials = [];
    for ii = 1:length(log_iv.Tone.tstart)

        bl = restrict(licks, log_iv.Tone.tstart(ii)-20, log_iv.Tone.tstart(ii));
        trials.base.licks(ii) = length(bl.t{1})./20;

        tl = restrict(licks, log_iv.Tone.tstart(ii), log_iv.Tone.tend(ii));
        trials.tone.licks(ii) = length(tl.t{1})./20;

        trl = restrict(licks, log_iv.Trace.tstart(ii), log_iv.Trace.tend(ii));
        trials.trace.licks(ii) = length(trl.t{1})./15;

        % same thing but movement
        bl = restrict(wheel_tsd, log_iv.Tone.tstart(ii)-20, log_iv.Tone.tstart(ii));
        trials.base.mov(ii) = sum(bl.data(3,:)) / length(bl.tvec);

        tl = restrict(wheel_tsd, log_iv.Tone.tstart(ii), log_iv.Tone.tend(ii));
        trials.tone.mov(ii) = sum(tl.data(3,:)) / length(tl.tvec);

        trl = restrict(wheel_tsd, log_iv.Trace.tstart(ii), log_iv.Trace.tend(ii));
        trials.trace.mov(ii) = sum(trl.data(3,:)) / length(trl.tvec);

    end


    % hold over

    parts = strsplit(s_list(iS).name, '_');
    data_out{iS}.name = [parts{3} parts{4}];
    data_out{iS}.sess = [parts{end-1} parts{end}(1:end-4)];

    data_out{iS}.licks = licks;
    data_out{iS}.wheel_tsd = wheel_tsd;
    data_out{iS}.log_iv = log_iv;
    data_out{iS}.trials = trials;

    sub_list{iS} = data_out{iS}.name;
    sess_list{iS} = data_out{iS}.sess;

    % data_wheel.(sess) = [data_wheel.(sess); data_out{iS}.wheel_tsd.data(3,:)]
end


ylabel('mean licks rate')

set(gca, 'XTick', 1:3, 'XTickLabel', {'Baseline [-20:0]', 'Tone', 'Trace'})
cond_1_idx = contains(sess_list, 'Cond1');
cond_2_idx = contains(sess_list, 'Cond2');
cond_3_idx = contains(sess_list, 'Cond3');
% for ii = 1:
% D1_mov = mean(data_out(cond1_idx).wheel_tsd.data(3,:))

%% plot each subject

c_ord = MS_linspecer(9);
r_ord = winter(6);

for iS =1:length(data_out)
    figure(iS)
    clf

    ax(1) = subplot(2,2,[1 3]);
    cla
    hold on
    for ii = 1:length(data_out{iS}.log_iv.Tone.tstart)

        if contains(data_out{iS}.sess, 'Rec')
            this_r = restrict(data_out{iS}.wheel_tsd, data_out{iS}.log_iv.Tone.tstart(ii)-20, data_out{iS}.log_iv.Trace.tend(ii)+20);
            this_l = restrict(data_out{iS}.licks, data_out{iS}.log_iv.Tone.tstart(ii)-20, data_out{iS}.log_iv.Trace .tend(ii)+20);
        else
            this_r = restrict(data_out{iS}.wheel_tsd, data_out{iS}.log_iv.Tone.tstart(ii)-20, data_out{iS}.log_iv.Puff.tend(ii)+20);

            % add the licks as a raster
            this_l = restrict(data_out{iS}.licks, data_out{iS}.log_iv.Tone.tstart(ii)-20, data_out{iS}.log_iv.Puff.tend(ii)+20);

        end
        plot(this_r.tvec - this_r.tvec(1)-20, this_r.data(2,:), 'color', c_ord(ii,:));

        plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
    end

    xline([0 20 35 36])
    rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
    rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
    rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(4,:), 'EdgeColor','none' )
    if contains(data_out{iS}.sess, 'Con')
        rectangle(gca, 'Position',[35 -.25 1 .25], 'FaceColor',r_ord(6,:), 'EdgeColor','none' )
    end
    rectangle(gca, 'Position',[36 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

    ylabel('Movement (a.u.)')
    xlim([-20 56])
    y_l = ylim;


    % bar plots summarizing movement and licks.
    ax(2) = subplot(3,2, 2);
    hold on
    hb = MS_bar_w_err(data_out{iS}.trials.tone.mov, data_out{iS}.trials.base.mov, c_ord(2:3,:), 1, [],[1 0]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_out{iS}.trials.tone.mov, data_out{iS}.trials.trace.mov, c_ord(1:2,:), 1,[] ,[1 2]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % ')
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Baseline [-20:0]', 'Tone', 'Trace'})

    % licks mean.
    ax(2) = subplot(2,2, 4);
    hold on
    hb = MS_bar_w_err(data_out{iS}.trials.tone.licks, data_out{iS}.trials.base.licks, c_ord(2:3,:), 1, 'ttest',[1 0]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_out{iS}.trials.tone.licks, data_out{iS}.trials.trace.licks, c_ord(1:2,:), 1, 'ttest',[1 2]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mean licks rate')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'Tone', 'Trace'})


    title([data_out{iS}.name ' ' data_out{iS}.sess])
end



%%
%% over Conditioning sessions

s_list = dir([ctrl_dir filesep 'arduino_log*Rec*.csv']);


% loop
data__rec_out = [];
for iS  = length(s_list):-1:1
    log_tab = readtable([ctrl_dir s_list(iS).name]);

    log_tab.time = log_tab.time - log_tab.time(1);  % zero to the first time point.
    log_tab.time = log_tab.time/ 1000; % convert from ms to seconds;

    % convert to evts

    labels = unique(log_tab.phase);

    phases = {'baseline', 'tone1', 'tone2'};

    for ii = 1:length(phases)
        log.(phases{ii}) = log_tab.time(contains(log_tab.event, '-') & contains(log_tab.phase, phases{ii}));
    end

    log.end = log_tab.time(contains(log_tab.event, 'Finish'));

    all_phases.t = log_tab.time(contains(log_tab.event, '-'));
    all_phases.data = log_tab.phase(contains(log_tab.event, '-'));

    % make intervals for each phase state
    log_iv =[];
    for ii = 1:length(phases)
        idx = find(ismember(all_phases.t, log.(phases{ii})));
        if max(idx) == length(all_phases.t)
            log_iv.(phases{ii}) = iv(all_phases.t(idx), [all_phases.t(idx(1:end-1)+1); log.end]);
        else
            log_iv.(phases{ii}) = iv(all_phases.t(idx), all_phases.t(idx+1));
        end
    end

    % get the encoder changes as a vector.

    wheel.tvec = log_tab.time(contains(log_tab.event, 'pos'));
    wheel.data = log_tab.encoderCount(contains(log_tab.event, 'pos'));

    % interp
    tvec_i = 0:.01:log.end;
    wheel.data = interp1(wheel.tvec, wheel.data, tvec_i);
    wheel.data(isnan(wheel.data)) = 0;
    wheel.tvec = tvec_i;

    % remove extreme values
    spd = diff([wheel.data wheel.data(end)]);

    spd(abs(spd) > 10 ) = NaN;

    spd = fillmissing(spd, 'nearest');

    wheel_tsd = tsd(wheel.tvec, [wheel.data; spd], {'encoder' 'spd'});

    % licks
    licks = []; licks.type = 'ts';
    licks.t{1} = log_tab.time(contains(log_tab.event, 'Lick'));
    licks.label{1} = 'licks';
    licks.cfg.history.mfun = {'ts'};
    licks.cfg.history.cfg = {[]};


    mov_bin = wheel_tsd.data(2,:) > .1;
    wheel_tsd.data(3,:) = mov_bin;
    wheel_tsd.label{3} = 'move binary';

    % count the licks per trial
    trials = [];
    for ii = 1:length(log_iv.tone1.tstart)

        bl = restrict(licks, log_iv.tone1.tstart(ii)-20, log_iv.tone1.tstart(ii));
        trials.base.licks(ii) = length(bl.t{1})./20;

        tl = restrict(licks, log_iv.tone1.tstart(ii), log_iv.tone1.tend(ii));
        trials.tone.licks(ii) = length(tl.t{1})./20;

        trl = restrict(licks, log_iv.tone1.tend(ii), log_iv.tone1.tend(ii)+20);
        trials.trace.licks(ii) = length(trl.t{1})./20;

        % same thing but movement
        bl = restrict(wheel_tsd, log_iv.tone1.tstart(ii)-20, log_iv.tone1.tstart(ii));
        trials.base.mov(ii) = (sum(bl.data(3,:)) / length(bl.tvec))*100;

        tl = restrict(wheel_tsd, log_iv.tone1.tstart(ii), log_iv.tone1.tend(ii));
        trials.tone.mov(ii) = (sum(tl.data(3,:)) / length(tl.tvec))*100;

        trl = restrict(wheel_tsd, log_iv.tone1.tend(ii), log_iv.tone1.tend(ii)+20);
        trials.trace.mov(ii) = (sum(trl.data(3,:)) / length(trl.tvec))*100;

    end

    for ii = 1:length(log_iv.tone2.tstart)

        bl = restrict(licks, log_iv.tone2.tstart(ii)-20, log_iv.tone2.tstart(ii));
        trials.base2.licks(ii) = length(bl.t{1})./20;

        tl = restrict(licks, log_iv.tone2.tstart(ii), log_iv.tone2.tend(ii));
        trials.tone2.licks(ii) = length(tl.t{1})./20;

        trl = restrict(licks, log_iv.tone2.tend(ii), log_iv.tone2.tend(ii)+20);
        trials.trace2.licks(ii) = length(trl.t{1})./20;

        % same thing but movement
        bl = restrict(wheel_tsd, log_iv.tone2.tstart(ii)-20, log_iv.tone2.tstart(ii));
        trials.base2.mov(ii) = (sum(bl.data(3,:)) / length(bl.tvec))*100;

        tl = restrict(wheel_tsd, log_iv.tone2.tstart(ii), log_iv.tone2.tend(ii));
        trials.tone2.mov(ii) = (sum(tl.data(3,:)) / length(tl.tvec))*100;

        trl = restrict(wheel_tsd, log_iv.tone2.tend(ii), log_iv.tone2.tend(ii)+20);
        trials.trace2.mov(ii) = (sum(trl.data(3,:)) / length(trl.tvec))*100;

    end


    % hold over
    parts = strsplit(s_list(iS).name, '_');
    data_rec_out{iS}.name = [parts{3} parts{4}];
    data_rec_out{iS}.sess = [parts{end-1} parts{end}(1:end-4)];

    data_rec_out{iS}.licks = licks;
    data_rec_out{iS}.wheel_tsd = wheel_tsd;
    data_rec_out{iS}.log_iv = log_iv;
    data_rec_out{iS}.trials = trials;

    sub_list{iS} = data_rec_out{iS}.name;
    sess_list{iS} = data_rec_out{iS}.sess;

    % data_wheel.(sess) = [data_wheel.(sess); data_out{iS}.wheel_tsd.data(3,:)]
end

%% load some blink tracking.

if ismac
    eye_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/Front_videos';
elseif ispc
    eye_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\Front_videos';
end
cd(eye_dir)

s_idx = 71; %for HF1 Rec

pos = MS_DLC2TSD_single('HF_1_Rec_d1DLC_Resnet50_HF_TFC_wide_cameraSep5shuffle1_snapshot_090_filtered.csv', 'HF_1_Rec_d1DLC_Resnet50_HF_TFC_wide_cameraSep5shuffle1_snapshot_090_filtered_p60_labeled.mp4', [1 1], 0);

pos = restrict(pos, pos.tvec(s_idx), pos.tvec(s_idx)+2550);

pos.tvec = pos.tvec - pos.tvec(1);

eye_d =[]; eye_v = [];

eye_d(1,:) = smoothdata(pos.data(1,:)', 'movmedian', 20);
eye_d(2,:) = smoothdata(pos.data(2,:)', 'movmedian', 20);

eye_v(1,:) = smoothdata(pos.data(5,:)', 'movmedian', 20);
eye_v(2,:) = smoothdata(pos.data(6,:)', 'movmedian', 20);

% d = pdist(pos.data(1:2,:)', pos.data(5:6,:)','euclidean');
d = [];
for ii = length(pos.data(1,:)):-1:1

    d(ii) = norm(eye_d(:,ii) - eye_v(:,ii));

end

d_z = zscore(d);

rm_idx =  abs(d_z) > 2.5 ;

d_z(rm_idx) = NaN;

d = fillmissing(d, 'nearest');

eye_tsd = tsd(pos.tvec, zscore(d), 'z eye dist');


for ii = 1:length(data_rec_out{1}.log_iv.tone1.tstart)

    bl = restrict(eye_tsd, data_rec_out{1}.log_iv.tone1.tstart(ii)-20, data_rec_out{1}.log_iv.tone1.tstart(ii));
    data_rec_out{1}.trials.base.eye(ii) = mean(bl.data);

    tl = restrict(eye_tsd, data_rec_out{1}.log_iv.tone1.tstart(ii), data_rec_out{1}.log_iv.tone1.tend(ii));
    data_rec_out{1}.trials.tone.eye(ii) = mean(tl.data);

    trl = restrict(eye_tsd, data_rec_out{1}.log_iv.tone1.tend(ii), data_rec_out{1}.log_iv.tone1.tend(ii)+20);
    data_rec_out{1}.trials.trace.eye(ii) = mean(trl.data);

end

for ii = 1:length(data_rec_out{1}.log_iv.tone2.tstart)

    bl = restrict(eye_tsd, data_rec_out{1}.log_iv.tone2.tstart(ii)-20, data_rec_out{1}.log_iv.tone2.tstart(ii));
    data_rec_out{1}.trials.base2.eye(ii) = mean(bl.data);

    tl = restrict(eye_tsd, data_rec_out{1}.log_iv.tone2.tstart(ii), data_rec_out{1}.log_iv.tone2.tend(ii));
    data_rec_out{1}.trials.tone2.eye(ii) = mean(tl.data);

    trl = restrict(eye_tsd, data_rec_out{1}.log_iv.tone2.tend(ii), data_rec_out{1}.log_iv.tone2.tend(ii)+20);
    data_rec_out{1}.trials.trace2.eye(ii) = mean(trl.data);

end


%% add some SWR data


for ii = 1:length(data_rec_out{1}.log_iv.tone1.tstart)

    bl = restrict(swr_r_iv, data_rec_out{1}.log_iv.tone1.tstart(ii)-20, data_rec_out{1}.log_iv.tone1.tstart(ii));
    data_rec_out{1}.trials.base.swr(ii) = length(bl.tstart)./20;

    tl = restrict(swr_r_iv, data_rec_out{1}.log_iv.tone1.tstart(ii), data_rec_out{1}.log_iv.tone1.tend(ii));
    data_rec_out{1}.trials.tone.swr(ii) = length(tl.tstart)./20;

    trl = restrict(swr_r_iv, data_rec_out{1}.log_iv.tone1.tend(ii), data_rec_out{1}.log_iv.tone1.tend(ii)+20);
    data_rec_out{1}.trials.trace.swr(ii) = length(trl.tstart)./20;

end

for ii = 1:length(data_rec_out{1}.log_iv.tone2.tstart)

    bl = restrict(swr_r_iv, data_rec_out{1}.log_iv.tone2.tstart(ii)-20, data_rec_out{1}.log_iv.tone2.tstart(ii));
    data_rec_out{1}.trials.base2.swr(ii) = length(bl.tstart)./20;

    tl = restrict(swr_r_iv, data_rec_out{1}.log_iv.tone2.tstart(ii), data_rec_out{1}.log_iv.tone2.tend(ii));
    data_rec_out{1}.trials.tone2.swr(ii) = length(tl.tstart)./20;

    trl = restrict(swr_r_iv, data_rec_out{1}.log_iv.tone2.tend(ii), data_rec_out{1}.log_iv.tone2.tend(ii)+20);
    data_rec_out{1}.trials.trace2.swr(ii) = length(trl.tstart)./20;

end
%% summar plots for Rec sessions

c_ord = winter(9);
r_ord = MS_linspecer(4);
s_ord = spring(9);

for iS =1%:length(data_rec_out)
    figure(iS)
    clf
    title(data_rec_out{iS}.name)

    ax(1) = subplot(3,6,1:3);
    cla
    hold on
    these_r = nan(length(data_rec_out{iS}.log_iv.tone1.tstart), 7501);
    these_eye = nan(length(data_rec_out{iS}.log_iv.tone1.tstart), 2250); 

    for ii = length(data_rec_out{iS}.log_iv.tone1.tstart):-1:1

        this_r = restrict(data_rec_out{iS}.wheel_tsd, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);
        this_l = restrict(swr_r_iv, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);
        this_eye = restrict(eye_tsd, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);

            these_eye(ii,1:length(this_eye.data)) = ((this_eye.data)); 

        these_r(ii,1:length(this_r.data)) = abs(this_r.data(2,:)); 

        % hdl(ii) = plot(this_r.tvec - this_r.tvec(1)-20, abs(this_r.data(2,:)), 'color', c_ord(ii,:), 'LineWidth',.25);

        plot([this_l.tstart- this_r.tvec(1)-20 this_l.tstart- this_r.tvec(1)-20]', [ones(size(this_l.tstart))*-(ii*.25)-.25 ones(size(this_l.tstart))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
    end

        % hdl = plot(this_r.tvec - this_r.tvec(1)-20, nanmean(these_r), 'color', c_ord(ii,:), 'LineWidth',1);
       hdl(1) = shadedErrorBar(this_r.tvec - this_r.tvec(1)-20, mean(these_r,"omitnan"), MS_SEM_vec(these_r), 'b');
            % plot(this_eye.tvec - this_eye.tvec(1)-20, (this_eye.data(1,:))+4, 'color', c_ord(ii,:), 'LineWidth',1);
            yyaxis right
       hdl(2) = shadedErrorBar(this_eye.tvec - this_eye.tvec(1)-20, mean(these_eye,"omitnan") , MS_SEM_vec(these_eye), 'g');
ylim([-8 2])

yyaxis left
    legend(gca, {'speed' 'eyelid diameter'}, 'box', 'off')

    xline([0 20 35])
    rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
    rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(1,:), 'EdgeColor','none' )
    rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
    rectangle(gca, 'Position',[35 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

    yline(2, '--')

    ylabel('speed (a.u)')
    xlim([-20 55])
ylim([-2 3])
title("CS+ | context 2")

    % bar plots summarizing movement and licks.
    subplot(3,6, 4);
    hold on
    hb = MS_bar_w_err(data_rec_out{iS}.trials.base.mov, data_rec_out{iS}.trials.tone.mov, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.mov, data_rec_out{iS}.trials.trace.mov, r_ord(1:2,:), 1,'ttest' ,[1 2]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % ')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
    ylim([0 70])

    subplot(3,6, 5);
    hold on
    hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.mov, data_rec_out{iS}.trials.tone2.mov, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_rec_out{iS}.trials.tone2.mov, data_rec_out{iS}.trials.trace2.mov, r_ord(1:2,:), 1,'ttest' ,[1 2]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % ')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', ' CS- trace'})
    ylim([0 70])


    subplot(3,6, 6);
    hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.mov, data_rec_out{iS}.trials.tone2.mov, [r_ord(1,:) ;r_ord(3,:)], 1,'ttest2' ,[0 1]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_rec_out{iS}.trials.trace.mov, data_rec_out{iS}.trials.trace2.mov, [r_ord(2,:) ;r_ord(4,:)], 1,'ttest2' ,[2 3]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('mov % ')
    set(gca, 'XTick', 0:4, 'XTickLabel', {'CS+', 'CS-',  'CS+ trace', 'CS- trace'})
    ylim([0 70])



    % swr mean.
    subplot(3,6, 10);
    cla
    hold on
    hb = MS_bar_w_err(data_rec_out{iS}.trials.base.swr, data_rec_out{iS}.trials.tone.swr, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.swr, data_rec_out{iS}.trials.trace.swr, r_ord(1:2,:), 1,'ttest' ,[1 2]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    %     hb = MS_bar_w_err(data_rec_out{iS}.trials.base.swr, data_rec_out{iS}.trials.trace.swr, r_ord(1:2,:), 1,'ttest' ,[0 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('swr rate ')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
    ylim([0 .5])


    subplot(3,6, 11);
        cla

    hold on
    hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.swr, data_rec_out{iS}.trials.tone2.swr, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_rec_out{iS}.trials.tone2.swr, data_rec_out{iS}.trials.trace2.swr, r_ord(1:2,:), 1,'ttest' ,[1 2]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    %     hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.swr, data_rec_out{iS}.trials.trace2.swr, [.7 .7 .7; r_ord(2,:)], 1,'ttest' ,[0 2]);
    % hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('swr rate ')
    set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
    ylim([0 .5])


    subplot(3,6, 12);
    hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.swr, data_rec_out{iS}.trials.tone2.swr, [r_ord(1,:) ;r_ord(3,:)], 1,'ttest2' ,[0 1]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    hb = MS_bar_w_err(data_rec_out{iS}.trials.trace.swr, data_rec_out{iS}.trials.trace2.swr, [r_ord(2,:) ;r_ord(4,:)], 1,'ttest2' ,[2 3]);
    hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

    ylabel('lick rate ')
    set(gca, 'XTick', 0:4, 'XTickLabel', {'CS+', 'CS-',  'C+ trace', 'CS- trace'})
    ylim([0 .5])

    % eyes
    if strcmp(data_rec_out{iS}.name, 'HF1')
        subplot(3,6, 16);
        hold on
        hb = MS_bar_w_err(data_rec_out{iS}.trials.base.eye, data_rec_out{iS}.trials.tone.eye, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
        hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

        hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.eye, data_rec_out{iS}.trials.trace.eye, r_ord(1:2,:), 1,'ttest' ,[1 2]);
        hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

        ylabel('z eyelid distance ')
        set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
        ylim([-1 2])

        subplot(3,6, 17);
        hold on
        hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.eye, data_rec_out{iS}.trials.tone2.eye, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
        hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

        hb = MS_bar_w_err(data_rec_out{iS}.trials.tone2.eye, data_rec_out{iS}.trials.trace2.eye, r_ord(1:2,:), 1,'ttest' ,[1 2]);
        hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

        ylabel('z eyelid distance ')
        set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
        ylim([-1 2])


        subplot(3,6, 18);
        hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.eye, data_rec_out{iS}.trials.tone2.eye, [r_ord(1,:) ;r_ord(3,:)], 1,'ttest2' ,[0 1]);
        hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

        hb = MS_bar_w_err(data_rec_out{iS}.trials.trace.eye, data_rec_out{iS}.trials.trace2.eye, [r_ord(2,:) ;r_ord(4,:)], 1,'ttest2' ,[2 3]);
        hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

        ylabel('z eyelid distance ')
        set(gca, 'XTick', 0:4, 'XTickLabel', {'CS+', 'CS-',  'C+ trace', 'CS- trace'})
        ylim([-1 2])
    end


    % tone 2

    ax(2) = subplot(3,6,7:9);
    cla
    hold on
    these_r = nan(length(data_rec_out{iS}.log_iv.tone2.tstart), 7501);
    these_eye = nan(length(data_rec_out{iS}.log_iv.tone2.tstart), 2250); 

    for ii = length(data_rec_out{iS}.log_iv.tone2.tstart):-1:1

        this_r = restrict(data_rec_out{iS}.wheel_tsd, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);
        this_l = restrict(swr_r_iv, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);
        this_eye = restrict(eye_tsd, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);

            these_eye(ii,1:length(this_eye.data)) = (this_eye.data); 

        these_r(ii,1:length(this_r.data)) = abs(this_r.data(2,:)); 

        % hdl(ii) = plot(this_r.tvec - this_r.tvec(1)-20, abs(this_r.data(2,:)), 'color', c_ord(ii,:), 'LineWidth',.25);

        plot([this_l.tstart- this_r.tvec(1)-20 this_l.tstart- this_r.tvec(1)-20]', [ones(size(this_l.tstart))*-(ii*.25)-.25 ones(size(this_l.tstart))*-(ii*.25)]', 'Color',s_ord(ii,:), 'LineWidth',2)
    end

        % hdl = plot(this_r.tvec - this_r.tvec(1)-20, nanmean(these_r), 'color', c_ord(ii,:), 'LineWidth',1);
       hdl(1) = shadedErrorBar(this_r.tvec - this_r.tvec(1)-20, mean(these_r,"omitnan"), MS_SEM_vec(these_r), 'r');
       hdl(1).mainLine.Color = s_ord(1,:);  hdl(1).mainLine.LineWidth = 1; 
       hdl(1).patch.FaceColor = s_ord(1,:)./[1 2 2]; 
       hdl(1).patch.FaceAlpha = .2; 
       hdl(1).mainLine.Color = s_ord(1,:)./[1 2 2]; 


yyaxis right
            % plot(this_eye.tvec - this_eye.tvec(1)-20, (this_eye.data(1,:))+4, 'color', c_ord(ii,:), 'LineWidth',1);
       hdl(2) = shadedErrorBar(this_eye.tvec - this_eye.tvec(1)-20, mean(these_eye,"omitnan") , MS_SEM_vec(these_eye),'r');
       hdl(2).mainLine.Color = s_ord(6,:);  hdl(2).mainLine.LineWidth = 1; 
       hdl(2).patch.FaceColor = s_ord(6,:)./[1 2 2]; 
       hdl(2).patch.FaceAlpha = .2; 
       hdl(2).mainLine.Color = s_ord(6,:)./[1 2 2]; 

ylim([-8 2])

yyaxis left
    legend(gca, {'speed' 'eyelid diameter'}, 'box', 'off')

    xline([0 20 35 36])
    rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
    rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(1,:), 'EdgeColor','none' )
    rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
    rectangle(gca, 'Position',[35 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )

    yline(2, '--')
    ylabel('speed (a.u)')
    xlim([-20 55])
    y_l = ylim;
    title("CS- | context 2")
end


%%

% 
% c_ord = winter(9);
% r_ord = MS_linspecer(4);
% s_ord = spring(9);
% 
% for iS =1%:length(data_rec_out)
%     figure(iS)
%     clf
%     title(data_rec_out{iS}.name)
% 
%     ax(1) = subplot(3,6,1:3);
%     cla
%     hold on
%     for ii = 1:length(data_rec_out{iS}.log_iv.tone1.tstart)
% 
%         this_r = restrict(data_rec_out{iS}.wheel_tsd, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);
%         this_l = restrict(data_rec_out{iS}.licks, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);
%         % if strcmp(data_rec_out{iS}.name, 'HF1')
%         %     this_eye = restrict(eye_tsd, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);
%         %     plot(this_eye.tvec - this_eye.tvec(1)-20, (this_eye.data(1,:)*.5)-4, 'color', c_ord(ii,:), 'LineWidth',1);
%         % else
%         %     this_eye = [];
%         % end
% 
%         hdl(ii) = plot(this_r.tvec - this_r.tvec(1)-20, abs(this_r.data(2,:)), 'color', c_ord(ii,:), 'LineWidth',1);
% 
%         plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
%     end
% 
%     legend(hdl([1 length(hdl)]), {'first trial' 'last trial'}, 'box', 'off')
% 
%     xline([0 20 35])
%     rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
%     rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(1,:), 'EdgeColor','none' )
%     rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
%     rectangle(gca, 'Position',[35 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
% 
%     ylabel('speed (a.u)')
%     xlim([-20 55])
%     y_l = ylim;
%     title("CS+ | context 2")
% 
%     % bar plots summarizing movement and licks.
%     subplot(3,6, 4);
%     hold on
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.base.mov, data_rec_out{iS}.trials.tone.mov, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.mov, data_rec_out{iS}.trials.trace.mov, r_ord(1:2,:), 1,'ttest' ,[1 2]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     ylabel('mov % ')
%     set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
%     ylim([0 70])
% 
%     subplot(3,6, 5);
%     hold on
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.mov, data_rec_out{iS}.trials.tone2.mov, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.tone2.mov, data_rec_out{iS}.trials.trace2.mov, r_ord(1:2,:), 1,'ttest' ,[1 2]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     ylabel('mov % ')
%     set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', ' CS- trace'})
%     ylim([0 70])
% 
% 
%     subplot(3,6, 6);
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.mov, data_rec_out{iS}.trials.tone2.mov, [r_ord(1,:) ;r_ord(3,:)], 1,'ttest2' ,[0 1]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.trace.mov, data_rec_out{iS}.trials.trace2.mov, [r_ord(2,:) ;r_ord(4,:)], 1,'ttest2' ,[2 3]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     ylabel('mov % ')
%     set(gca, 'XTick', 0:4, 'XTickLabel', {'CS+', 'CS-',  'CS+ trace', 'CS- trace'})
%     ylim([0 70])
% 
% 
% 
%     % licks mean.
%     subplot(3,6, 10);
%     hold on
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.base.licks, data_rec_out{iS}.trials.tone.licks, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.licks, data_rec_out{iS}.trials.trace.licks, r_ord(1:2,:), 1,'ttest' ,[1 2]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     ylabel('lick rate ')
%     set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
%     ylim([0 6])
% 
% 
%     subplot(3,6, 11);
%     hold on
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.licks, data_rec_out{iS}.trials.tone2.licks, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.tone2.licks, data_rec_out{iS}.trials.trace2.licks, r_ord(1:2,:), 1,'ttest' ,[1 2]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     ylabel('lick rate ')
%     set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
%     ylim([0 6])
% 
% 
%     subplot(3,6, 12);
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.licks, data_rec_out{iS}.trials.tone2.licks, [r_ord(1,:) ;r_ord(3,:)], 1,'ttest2' ,[0 1]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     hb = MS_bar_w_err(data_rec_out{iS}.trials.trace.licks, data_rec_out{iS}.trials.trace2.licks, [r_ord(2,:) ;r_ord(4,:)], 1,'ttest2' ,[2 3]);
%     hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%     ylabel('lick rate ')
%     set(gca, 'XTick', 0:4, 'XTickLabel', {'CS+', 'CS-',  'C+ trace', 'CS- trace'})
%     ylim([0 6])
% 
%     % eyes
%     if strcmp(data_rec_out{iS}.name, 'HF1')
%         subplot(3,6, 16);
%         hold on
%         hb = MS_bar_w_err(data_rec_out{iS}.trials.base.eye, data_rec_out{iS}.trials.tone.eye, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
%         hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%         hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.eye, data_rec_out{iS}.trials.trace.eye, r_ord(1:2,:), 1,'ttest' ,[1 2]);
%         hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%         ylabel('z eyelid distance ')
%         set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS+', 'CS+ trace'})
%         ylim([-1 2])
% 
%         subplot(3,6, 17);
%         hold on
%         hb = MS_bar_w_err(data_rec_out{iS}.trials.base2.eye, data_rec_out{iS}.trials.tone2.eye, [.7 .7 .7; r_ord(1,:)], 1, 'ttest' ,[0 1]);
%         hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%         hb = MS_bar_w_err(data_rec_out{iS}.trials.tone2.eye, data_rec_out{iS}.trials.trace2.eye, r_ord(1:2,:), 1,'ttest' ,[1 2]);
%         hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%         ylabel('z eyelid distance ')
%         set(gca, 'XTick', 0:2, 'XTickLabel', {'Baseline [-20:0]', 'CS-', 'CS- trace'})
%         ylim([-1 2])
% 
% 
%         subplot(3,6, 18);
%         hb = MS_bar_w_err(data_rec_out{iS}.trials.tone.eye, data_rec_out{iS}.trials.tone2.eye, [r_ord(1,:) ;r_ord(3,:)], 1,'ttest2' ,[0 1]);
%         hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%         hb = MS_bar_w_err(data_rec_out{iS}.trials.trace.eye, data_rec_out{iS}.trials.trace2.eye, [r_ord(2,:) ;r_ord(4,:)], 1,'ttest2' ,[2 3]);
%         hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
%         ylabel('z eyelid distance ')
%         set(gca, 'XTick', 0:4, 'XTickLabel', {'CS+', 'CS-',  'C+ trace', 'CS- trace'})
%         ylim([-1 2])
%     end
% 
% 
%     % tone 2
% 
%     ax(2) = subplot(3,6,7:9);
%     cla
%     hold on
%     for ii = 1:length(data_rec_out{iS}.log_iv.tone2.tstart)
% 
%         this_r = restrict(data_rec_out{iS}.wheel_tsd, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);
%         this_l = restrict(data_rec_out{iS}.licks, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);
% 
%         % if strcmp(data_rec_out{iS}.name, 'HF1')
%         %     this_eye = restrict(eye_tsd, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);
%         %     plot(this_eye.tvec - this_eye.tvec(1)-20, (this_eye.data(1,:)*.5)-4, 'color', c_ord(ii,:), 'LineWidth',1);
%         % else
%         %     this_eye = [];
%         % end
% 
%         hdl(ii) =  plot(this_r.tvec - this_r.tvec(1)-20, abs(this_r.data(2,:)), 'color', c_ord(ii,:), 'LineWidth',1);
% 
%         plot([this_l.t{1}- this_r.tvec(1)-20 this_l.t{1}- this_r.tvec(1)-20]', [ones(size(this_l.t{1}))*-(ii*.25)-.25 ones(size(this_l.t{1}))*-(ii*.25)]', 'Color',c_ord(ii,:), 'LineWidth',2)
%     end
% 
%     legend(hdl([1 length(hdl)]), {'first trial' 'last trial'}, 'box', 'off')
% 
%     xline([0 20 35 36])
%     rectangle(gca, 'Position',[-20 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
%     rectangle(gca, 'Position',[0 -.25 20 .25], 'FaceColor',r_ord(1,:), 'EdgeColor','none' )
%     rectangle(gca, 'Position',[20 -.25 15 .25], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
%     rectangle(gca, 'Position',[35 -.25 20 .25], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
% 
%     ylabel('speed (a.u)')
%     xlim([-20 55])
%     y_l = ylim;
%     title("CS- | context 2")
% 
% 
% 
%     ax(2) = subplot(3,6,13:15);
%     cla
%     hold on
%     for ii = 1:length(data_rec_out{iS}.log_iv.tone2.tstart)
% 
%         if strcmp(data_rec_out{iS}.name, 'HF1')
%             yyaxis left
%             this_eye = restrict(eye_tsd, data_rec_out{iS}.log_iv.tone1.tstart(ii)-20, data_rec_out{iS}.log_iv.tone1.tend(ii)+35);
%             plot(this_eye.tvec - this_eye.tvec(1)-20, (this_eye.data(1,:)), '-','color', c_ord(ii,:), 'LineWidth',1);
% 
%             yyaxis right
%             this_eye = restrict(eye_tsd, data_rec_out{iS}.log_iv.tone2.tstart(ii)-20, data_rec_out{iS}.log_iv.tone2.tend(ii)+35);
%             plot(this_eye.tvec - this_eye.tvec(1)-20, (this_eye.data(1,:)-6),'-',  'color', s_ord(ii,:), 'LineWidth',1);
%         else
%             this_eye = [];
%         end
% 
%     end
%     yyaxis left
%     legend(hdl([1 length(hdl)]), {'first trial' 'last trial'}, 'box', 'off')
% 
%     xline([0 20 35 36])
%     rectangle(gca, 'Position',[-20 -3 20 1], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
%     rectangle(gca, 'Position',[0 -3 20 1], 'FaceColor',r_ord(1,:), 'EdgeColor','none' )
%     rectangle(gca, 'Position',[20 -3 15 1], 'FaceColor',r_ord(2,:), 'EdgeColor','none' )
%     rectangle(gca, 'Position',[35 -3 20 1], 'FaceColor',[.7 .7 .7], 'EdgeColor','none' )
% 
% 
%     ylabel('eye dist (zscore)')
%     xlim([-20 55])
%     yyaxis right
% 
%     ylim([-10 6]);
%     set(gca,'ytick', -8:2:8)
%     set(gca, "YTickLabel", get(gca, 'ytick') +6)
% 
%     yyaxis left
%     ylim([-10 6]);
% 
%     text(-15, 5, 'CS+', 'HorizontalAlignment','left', 'FontSize',14)
% 
%     text(-15, -9, 'CS-', 'HorizontalAlignment','left', 'FontSize',14)
% 
%     % title("CS- | context 2")
% end
