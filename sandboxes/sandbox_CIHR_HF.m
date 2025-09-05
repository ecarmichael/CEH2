%% sandbox_CIHR_HR_TFC

if ispc
    kilo_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-03_16-34-57_TFC_ENC\Record Node 113\experiment1\recording1\continuous\Intan_RHD_USB-100.Rhythm Data';
    OE_evts_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-03_16-34-57_TFC_ENC\Record Node 113\experiment1\recording1\events\Intan_RHD_USB-100.Rhythm Data\TTL';
    ctrl_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\conditioning 1/';
    eye_dir = [];
elseif ismac
    kilo_dir = [];
    ctrl_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/conditioning 1/';
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
wheel.data(3,:) = mov_bin; 
wheel.label{3} = 'move binary';

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

hb = MS_bar_w_err(trials.tone.licks, trials.base.licks, c_ord(2:3,:), 1, 'ttest',[1 0]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

hb = MS_bar_w_err(trials.tone.licks, trials.trace.licks, c_ord(1:2,:), 1, 'ttest',[1 2]); 

hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

ylabel('mean licks rate')

set(gca, 'XTick', 1:3, 'XTickLabel', {'Baseline [-20:0]', 'Tone', 'Trace'})