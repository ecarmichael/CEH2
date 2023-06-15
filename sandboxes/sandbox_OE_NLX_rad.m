%% sandbox Radial Ephys

% raw_data_dir = 'J:\M29_2023-03-31_12-51-20_Rad5\Record Node 113\experiment1\recording2\continuous\Acquisition_Board-100.Rhythm Data';
% 
% spike_dir = 'J:\M29_2023_04_01_Rad5_kilo'%'C:\Users\ecarm\Desktop\M29_2023_04_01_Rad5_kilo';
% 
% nlx_data_dir = 'D:\M29_2023-03-31_Rad5';

raw_data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Radial_ephys\M29\M29_2023-04-14_12-13-30_Rad10\Record Node 113\experiment1\recording2\continuous\Acquisition_Board-100.Rhythm Data';
oe_evt_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Radial_ephys\M29\M29_2023-04-14_12-13-30_Rad10\Record Node 113\experiment1\recording2\events\Acquisition_Board-100.Rhythm Data\TTL'; 
spike_dir = 'C:\Users\ecarm\Desktop\M29_2023_04_14_Rad10_kilo\Kilo';%'C:\Users\ecarm\Desktop\M29_2023_04_01_Rad5_kilo';

nlx_data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Radial_ephys\M29\M29_2023-04-14_12-13-30_Rad10\M29_2023-04-14_Rad10';

csc_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Radial_ephys\M29\M29_2023-04-14_12-13-30_Rad10\Record Node 123'; 
c_ord = linspecer(7); 

%% load the events
cd(oe_evt_dir)
oe_evt.states = readNPY('states.npy');
oe_evt.t = readNPY('timestamps.npy');

oe_t_start = oe_evt.t(oe_evt.states == 3); 
s_idx = find((diff(oe_t_start) > 5)); 
oe_t_start = oe_t_start([1; s_idx(1:end)+1]);

% oe_t_start_s = oe_evt.states(oe_evt.states == 4); 
% oe_t_start_s = oe_t_start_s([1; s_idx(1:end)+1]);

oe_t_end = oe_evt.t(oe_evt.states == 4); 
e_idx = find((diff(oe_t_end) > 5)); 
oe_t_end = oe_t_end([1; e_idx(1:end)+1]);
% oe_t_end = [oe_t_end(1) - 60; oe_t_end]; 

% oe_t_end_s = oe_evt.states(oe_evt.states == 3); 
% oe_t_end_s = oe_t_end_s([1; e_idx(1:end)+1]);

% check everything
% figure(1919)
% cla
% hold on
% plot(oe_t_start, ones(size(oe_t_start)), 'dr')
% plot(oe_t_end, ones(size(oe_t_end)), 'sb')
% legend({'start', 'end'})

%% get some tvec info
cd(raw_data_dir)
tvec = readNPY('timestamps.npy');

% tvec = oe_evt.t(1) - tvec(1); 
%% load the tracking

cd(nlx_data_dir);

evt = LoadEvents([]);


start_rec= evt.t{ismember(evt.label, 'Starting Recording')};
end_rec = evt.t{ismember(evt.label, 'Stopping Recording')};

sleep_times = evt.t{ismember(evt.label, 'sleep')};
% sleep_times = sleep_times(1);

if sum(ismember(evt.label, 'end sleep')) == 0
    if length(sleep_times) == 2 && diff(sleep_times) > 3.5*60*60
        fprintf('Two sleep times found. Duration: %5.2f sec (%5.2f hrs) \n', diff(sleep_times)/60, diff(sleep_times)/60/60)
        
        end_sleep_time = sleep_times(end);
        start_sleep_time = sleep_times(1);
        %     end_sleep_times =   sleep_times + 4*60*60;
    end
else
    end_sleep_times = evt.t{ismember(evt.label, 'end sleep')};
end

if sum(ismember(evt.label, 'rem')) == 0
    REM_sleep_times = [];
else
    REM_sleep_times = evt.t{ismember(evt.label, 'rem')};
end

% trial events

% end_sleep_times = evt.t{ismember(evt.label, 'end sleep')};
% REM_sleep_times = evt.t{ismember(evt.label, 'rem')};

t_start = [start_rec evt.t{ismember(evt.label, 't')}];
t_end = [evt.t{ismember(evt.label, 'te')}];

% normalize to start rec
end_rec = end_rec - start_rec;

start_sleep_time = start_sleep_time - start_rec;
end_sleep_time = end_sleep_time - start_rec;

REM_sleep_times = REM_sleep_times - start_rec;

t_start = t_start - start_rec;
t_end = t_end - start_rec;

start_rec_og = start_rec;
start_rec = start_rec - start_rec;


% fix overlapping t and te
over_lap_idx = ismember(t_start, t_end);
t_start(over_lap_idx) = [];

t_start_encode = t_start(t_start < start_sleep_time);
t_start_recall = t_start(t_start > end_sleep_time);

t_end_encode = t_end(t_end < start_sleep_time);
t_end_recall = t_end(t_end > end_sleep_time);

%% align to first trial end. 
% nlx_t_end_1 = t_end(1); 

% oe_t_end_1 = oe_t_end(1) ;

figure(10101)
cla
hold on
scatter(oe_t_start, ones(size(oe_t_start)),50, 'sb', 'filled')
scatter(oe_t_end, ones(size(oe_t_end)),50, 'db', 'filled')
scatter(t_start, ones(size(t_start))+1,50, 'sr', 'filled')
scatter(t_end, ones(size(t_end))+1,50, 'dr', 'filled')
ylim([0 2.5])
% plot(t_end, ones(size(t_end)), 'dr')

offset = mean(oe_t_end - t_end'); 
fprintf('estimated offset: %0.2fs. Correcting... \n', offset)

scatter(oe_t_start - offset, ones(size(oe_t_start))+.5,50, c_ord(3,:), 's', 'filled')
scatter(oe_t_end- offset, ones(size(oe_t_end))+.5,50, c_ord(3,:), 'd', 'filled')
% scatter(t_start, ones(size(t_start))+0.001, 'sr', 'filled')
% scatter(t_end, ones(size(t_end))+0.001, 'dr', 'filled')
legend({ 'oe  start', 'oe end','nlx start', 'nlx end', 'oe correct start', 'oe corrected end'})

set(gca, 'ytick', [ 1 1.5 2], 'yticklabels', {'OE org', 'OE Corrected', 'NLX'})


%% load the position and trim it to the trial and sleep periods. 
pos_cfg.convFact = [6.4 6.4];
pos = LoadPos(pos_cfg);



% apply light smoothing 
pos.data(1,:) = smoothdata(pos.data(1,:),'gaussian', round(1/mode(diff(pos.tvec)))/2); 
pos.data(2,:) = smoothdata(pos.data(2,:),'gaussian', round(1/mode(diff(pos.tvec)))/2); 

linspeed = getLinSpd([],pos); % linear speed
 
% Threshold speed
cfg = []; cfg.method = 'raw'; cfg.operation = 'range'; cfg.threshold = [4 20]; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg,linspeed); % only keep intervals with speed above thresh
 
% pos_fast = restrict(pos, iv_fast);
pos_fast = pos; 

pos_fast.tvec = pos_fast.tvec - pos.tvec(1);
pos.tvec = pos.tvec - pos.tvec(1);

pos_encode = restrict(pos, start_rec, start_sleep_time);
pos_recall = restrict(pos, end_sleep_time, pos.tvec(end));
pos_sleep = restrict(pos, start_sleep_time, end_sleep_time);

pos_encode_trials = restrict(pos_fast, t_start_encode+60, t_end_encode);
pos_encode_iti = restrict(pos_fast, t_start_encode, t_start_encode+60);


pos_recall_trials = restrict(pos_fast, t_start_recall(2:end)+60, t_end_recall(2:end));
pos_recall_iti = restrict(pos_fast, t_start_recall, t_start_recall+60);


fprintf('Encoding dur: %.0fmins \nSleep dur: %.0fmins \nRecall dur: %.0fmins\n', (pos_encode.tvec(end) - pos_encode.tvec(1))/60,  (pos_sleep.tvec(end) - pos_sleep.tvec(1))/60, (pos_recall.tvec(end) - pos_recall.tvec(1))/60)

clear pos

% test plot
figure(101)
clf
subplot(3,3,1)
hold on
plot(pos_encode_trials.data(1,:), pos_encode_trials.data(2,:), '.', 'color', c_ord(2,:));
plot(pos_encode_iti.data(1,:), pos_encode_iti.data(2,:), '.', 'color', c_ord(1,:));
legend({'Trials'; 'ITI'})
title('Ecoding phase')


%plot(pos_encode_trials.data(1,:), pos_encode_trials.data(2,:), '.', 'color', c_ord(1,:));
%plot(pos_encode_iti.data(1,:), pos_encode_iti.data(2,:), '.', 'color', c_ord(2,:));
%legend({'Trials'; 'ITI'})wake_t = [0 1660 11177 13440 17000 csc.tvec(end)]; 

%title('Ecoding phase')

subplot(3,3,2)
plot(pos_sleep.data(1,:), pos_sleep.data(2,:), '.k');
title('Sleep')

subplot(3,3,3)
hold on
plot(pos_recall_trials.data(1,:), pos_recall_trials.data(2,:), '.', 'color', c_ord(2,:));
plot(pos_recall_iti.data(1,:), pos_recall_iti.data(2,:), '.', 'color', c_ord(1,:));
title('Recall phase')


%% time position movie for debugginh

% figure(1010)
% cla

% for ii = 1:length(pos_recall.



%% load some spikes and convert them to the TS format.
cd(spike_dir)

S = OE_phy2TS;

% convert S.t to time using timestamps.
for ii = 1:length(S.t)
    this_t = tvec(S.t{ii}) - offset;
    S.t{ii} = this_t(this_t >0); 

end

%restrict to encoding/sleep/recall
S_encode = restrict(S, start_rec, start_sleep_time);
S_recall = restrict(S, end_sleep_time, end_rec);
S_sleep = restrict(S, start_sleep_time, end_sleep_time);

S_encode_trials = restrict(S, t_start_encode+60, t_end_encode);
S_encode_iti = restrict(S, t_start_encode, t_start_encode+60);

S_recall_trials = restrict(S, t_start_recall(2:end)+60, t_end_recall(2:end));
S_recall_iti = restrict(S, t_start_recall, t_start_recall+60);

%% plot the spike rasters

figure(101)
subplot(3,3,7:9)
plot(S);
% for ii = 1:length(REM_sleep_times)
% xline([REM_sleep_times(ii)], 'color', c_ord(2,:))
% end
xline([start_sleep_time], 'color', c_ord(3,:), 'linewidth', 5)
xline([end_sleep_time], 'color', c_ord(3,:), 'linewidth', 5)




subplot(3,3,4)
plot(S_encode);
for ii = 1:length(t_start_encode)
xline([t_start_encode(ii)], 'color', c_ord(2,:), 'linewidth', 5)
end
for ii = 1:length(t_end_encode)
xline([t_end_encode(ii)], 'color', c_ord(1,:), 'linewidth', 5)
end
% vline(t_start_encode, 'g')
% vline(t_end_encode, 'm')

subplot(3,3,5)
plot(S_sleep);
% vline(t_start_encode- start_rec, 'g')
% vline(t_end_encode- start_rec, 'm')

subplot(3,3,6)
plot(S_recall);
for ii = 1:length(t_start_recall)
xline([t_start_recall(ii)], 'color', c_ord(2,:), 'linewidth', 5)
end
for ii = 1:length(t_end_recall)
xline([t_end_recall(ii)], 'color', c_ord(1,:), 'linewidth', 5)
end
%% get the occupancy and prep for rate maps


% set up bins
SET_xmin = 0; SET_ymin = 20; % set up bins
SET_xmax = 80; SET_ymax = 110;
SET_xBinSz = 2.5; SET_yBinSz = 2.5;


 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([1 1],1); % 2d gaussian in bins

% compute occupancy encode
occ_hist_encode = histcn(pos_encode_trials.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist_encode = conv2(occ_hist_encode,kernel,'same');
 
no_occ_idx_encode = find(occ_hist_encode < 15); % NaN out bins never visited
occ_hist_encode(no_occ_idx_encode) = NaN;
 
occ_hist_encode = occ_hist_encode .* (1/30); % convert samples to seconds using video frame rate (30 Hz)

% compute occupancy recall
occ_hist_recall = histcn(pos_recall_trials.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist_recall = conv2(occ_hist_recall,kernel,'same');

no_occ_idx_recall = find(occ_hist_recall < 7.5); % NaN out bins visited for less than 1 s
occ_hist_recall(no_occ_idx_recall) = NaN;
 
occ_hist_recall = occ_hist_recall .* (1/30); % convert samples to seconds using video frame rate (30 Hz)



%% plot some simple maps with spikes
m = 4;
n = 4;
s1_idx = 1:2:n*m; 
s2_idx = 2:2:n*m; 
figure(102)
clf
ip = 0; 
for ii = 1:length(S.t)
    spk_x = interp1(pos_encode_trials.tvec,pos_encode_trials.data(1,:),S_encode_trials.t{ii},'linear');
    spk_y = interp1(pos_encode_trials.tvec,pos_encode_trials.data(2,:),S_encode_trials.t{ii},'linear');
    if ip >= (n*m)/2
        figure(102+ii)
        ip = 1;
    else
        ip = ip+1;
    end
    
    disp(ip)
    subplot(m, n, s1_idx(ip))
    plot(pos_encode_trials.data(1,:), pos_encode_trials.data(2,:), '.k');
    title(S.label{ii})
    hold on
    plot(spk_x, spk_y, '.r');axis off;
    
    % compute the rate map
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx_encode) = NaN;
    tc = spk_hist./occ_hist_encode;

    subplot(m, n, s2_idx(ip))
    pcolor(tc'); shading flat; axis off; %colorbar('Location', 'northoutside')
   
end


%% recall
%% plot some simple maps with spikes
m = 4;
n = 4;
s1_idx = 1:2:n*m; 
s2_idx = 2:2:n*m; 
figure(102)
clf
ip = 0; 
for ii = 1:length(S.t)
    spk_x = interp1(pos_recall_trials.tvec,pos_recall_trials.data(1,:),S_recall_trials.t{ii},'linear');
    spk_y = interp1(pos_recall_trials.tvec,pos_recall_trials.data(2,:),S_recall_trials.t{ii},'linear');
    if ip >= (n*m)/2
        figure(2002+ii)

        ip = 1;
    else
        ip = ip+1;
    end
    
    disp(ip)
    subplot(m, n, s1_idx(ip))
    plot(pos_recall_trials.data(1,:), pos_recall_trials.data(2,:), '.k');
    
    hold on
    plot(spk_x, spk_y, '.r');axis off;
        title(S.label{ii})

    % compute the rate map
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx_recall) = NaN;
    tc = spk_hist./occ_hist_recall;

    subplot(m, n, s2_idx(ip))
    pcolor(tc'); shading flat; axis off; %colorbar('Location', 'northoutside')
   
end


%% pick certain cells

cell_id = 104; 
ii = contains(S.label, num2str(cell_id))
% ii = 48; 
% ii = 33


figure(999)
subplot(2,2,1)
 spk_x = interp1(pos_encode_trials.tvec,pos_encode_trials.data(1,:),S_encode_trials.t{ii},'linear');
    spk_y = interp1(pos_encode_trials.tvec,pos_encode_trials.data(2,:),S_encode_trials.t{ii},'linear');
  plot(pos_encode_trials.data(1,:), pos_encode_trials.data(2,:), '.k');
    
    hold on
    plot(spk_x, spk_y, '.', 'color', c_ord(1,:));axis off;
        title([num2str(S.label{ii}) ' Encoding'])
        
        
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx_recall) = NaN;
    tc = spk_hist./occ_hist_recall;

    subplot(2,2,2)
    pcolor(tc'); shading flat; axis off;
        
        subplot(2,2,3)
         spk_x = interp1(pos_recall_trials.tvec,pos_recall_trials.data(1,:),S_recall_trials.t{ii},'linear');
    spk_y = interp1(pos_recall_trials.tvec,pos_recall_trials.data(2,:),S_recall_trials.t{ii},'linear');
  plot(pos_recall_trials.data(1,:), pos_recall_trials.data(2,:), '.k');
    
    hold on
    plot(spk_x, spk_y, '.', 'color', c_ord(2,:));axis off;
        title([num2str(S.label{ii}) ' Recall'])
        
        
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx_recall) = NaN;
    tc = spk_hist./occ_hist_recall;

    subplot(2,2,4)
    pcolor(tc'); shading flat; axis off
    
    
    %% try some SWR dtection
    cd(csc_dir)
    [data, tvec_csc, hdr] = load_open_ephys_data('100_RhythmData-B_CH9.continuous');

cfg.decimateByFactor = 30000 / 1000;
data = decimate(data,cfg.decimateByFactor);
tvec_csc = tvec_csc(1:cfg.decimateByFactor:end);

csc = tsd(tvec_csc - offset, data', {'CSC7.ncs'});
csc.cfg.hdr{1}.SamplingFrequency = 1000;
clear data tvec_csc


csc_sleep = restrict(csc,0,start_sleep_time)%, csc.tvec(end))

%% emg
[emg_data, emg_tvec, emg_hdr] = load_open_ephys_data('100_RhythmData-B_CH52.continuous');

emg_data = decimate(emg_data,cfg.decimateByFactor);
emg_tvec = emg_tvec(1:cfg.decimateByFactor:end);

emg = tsd(emg_tvec, emg_data', 'EMG');
emg.cfg.hdr{1}.SamplingFrequency = 1000;
clear emg_data emg_tvec

% get the theta delta ratio
cfg_emg = [];
cfg_emg.threshold = 0;
cfg_emg.f = [15 30]; 
% cfg_emg.type = 'cheby1';
% cfg_emg.display_filter = 1
emg_f = FilterLFP(cfg_emg,csc);

% emg_sleep = restrict(emg_f, end_sleep_time, csc.tvec(end))
emg_sleep = restrict(emg_f,0,start_sleep_time)%, csc.tvec(end))


%% movement
emg_rms = emg_sleep; 

emg_rms.data = sqrt(movmedian(emg_sleep.data.^2, emg_sleep.cfg.hdr{1}.SamplingFrequency*10));      % RMS Value Over ‘WinLen’ Samples

% emg_rms_z = zscore(emg_rms.data); 

emg_rms_mm = movmedian(emg_rms.data,(csc.cfg.hdr{1}.SamplingFrequency)*1);

move_thresh = prctile(emg_rms_mm, 10);
% REM_thresh = prctile(emg_rms_mm, 15); 

move_idx = emg_rms_mm >move_thresh; 
% move_ iv
% cfg = []; cfg.method = 'zscore'; cfg.operation = '>';cfg.dcn = cfg.operation ; cfg.threshold = .5; % speed limit in cm/sec
% iv_move = TSDtoIV(cfg,emg_rms); % only keep intervals with speed above thresh



%%
cfg_swr = [];
cfg_swr.check = 0; % plot checks.
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [125 200]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter = 0; % use this to see the fvtool

% smoothing
cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;

% detection
% cfg_swr.artif_det.method = 'zscore';
% cfg_swr.artif_det.threshold = 8;
% cfg_swr.artif_det.dcn = '>';
% cfg_swr.artif_det.rm_len = .2;
cfg_swr.threshold = 1.5;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
cfg_swr.merge_thr = 0.02; %merge events that are within 20ms of each other.

% 
% cfg_swr.nan_idx = nan_idx; % where are any nans, say from excluding artifacts, other events...

% restrictions
cfg_swr.max_len = [];
cfg_swr.max_len.operation = '<';
cfg_swr.max_len.threshold = .1;

%                 cfg_swr.min_len = [];
%                 cfg_swr.min_len.operation = '<';
%                 cfg_swr.min_len.threshold = .2;
cfg_swr.nCycles = 5; % number of cycles
cfg_swr.nCycles_operation = '>='; % number of cycles

% variaence
% cfg_swr.var = [];
% cfg_swr.var.operation = '<';
% cfg_swr.var.threshold = 1;
cfg_swr.nan_idx = move_idx;

SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc_sleep);

% remove periods of movement
% SWR_evts = DifferenceIV([], SWR_evts, mov_evts);

%
%         cfg_max_len = [];
%         cfg_max_len.operation = '>';
%         cfg_max_len.threshold = 5;
%          SWR_evts = SelectIV(cfg_max_len,SWR_evts,'nCycles');

%% keep only SWRs during imobility

% swr_keep = DifferenceIV([],SWR_evts, iv_move ); 

%manual event curration
for ii = 1:length(SWR_evts.tstart)
   SWR_evts.usr.id(ii) = ii; 
    
end

cfg_plot = [];
cfg_plot.display = 'iv'; %'iv';
cfg_plot.title = 'id';
PlotTSDfromIV(cfg_plot, swr_keep, csc_sleep)

%%
% keep_idx = [1 29 30 32 37 38 39 61 72 73 74 75 81 82 83 96 97 98 99 100 101 102 103 104 105 106 113 114 116 117 119 123 137 138 139 140 141 142 143 145 146 160 161 163 164 173 174 175 176 177 178 179 180 181 186 ];

keep_idx = [9 10 11 12 13 14 17 18 19 20 21 22 23 24 25 26 27];


rm_idx = ones(1,length(swr_keep.tstart)); 
rm_idx(keep_idx) = 0; 
rm_idx = logical(rm_idx); 
% for ii = length(swr_keep.tstart):-1:1
%     if ismember(ii, keep_idx)
        swr_keep.tstart(rm_idx) = [];
        swr_keep.tend(rm_idx) = [];
        f_names = fieldnames(swr_keep.usr);
        for iF = 1:length(f_names)
            swr_keep.usr.(f_names{iF})(rm_idx) = [];
        end
%     end
% end

subplot(2,1,1)
cla
% check quality. 
cfg_plot = [];
% cfg_plot.display = 'iv'; %'iv';
cfg_plot.title = 'id';
PlotTSDfromIV(cfg_plot, swr_keep, csc_sleep)
    


    