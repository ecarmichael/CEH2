%% sandbox Radial Ephys

raw_data_dir = 'J:\M29_2023-03-31_12-51-20_Rad5\Record Node 113\experiment1\recording2\continuous\Acquisition_Board-100.Rhythm Data';

spike_dir = 'J:\M29_2023_04_01_Rad5_kilo'%'C:\Users\ecarm\Desktop\M29_2023_04_01_Rad5_kilo';

nlx_data_dir = 'D:\M29_2023-03-31_Rad5';

c_ord = linspecer(3); 
%% load the tracking

cd(nlx_data_dir);

evt = LoadEvents([]);


start_rec= evt.t{ismember(evt.label, 'Starting Recording')};
end_rec = evt.t{ismember(evt.label, 'Stopping Recording')};

sleep_times = evt.t{ismember(evt.label, 'sleep')};
end_sleep_times = evt.t{ismember(evt.label, 'end sleep')};
REM_sleep_times = evt.t{ismember(evt.label, 'rem')};

t_start = [start_rec evt.t{ismember(evt.label, 't')}];
t_end = [evt.t{ismember(evt.label, 'te')}];

% normalize to start rec
end_rec = end_rec - start_rec;

sleep_times = sleep_times - start_rec;
end_sleep_times = end_sleep_times - start_rec;
REM_sleep_times = REM_sleep_times - start_rec;

t_start = t_start - start_rec;
t_end = t_end - start_rec;

start_rec = start_rec - start_rec;


% fix overlapping t and te
over_lap_idx = ismember(t_start, t_end);
t_start(over_lap_idx) = [];

t_start_encode = t_start(t_start < sleep_times);
t_start_recall = t_start(t_start > end_sleep_times);

t_end_encode = t_end(t_end < sleep_times);
t_end_recall = t_end(t_end > end_sleep_times);

pos = LoadPos([]);
pos.tvec = pos.tvec - pos.tvec(1);

pos_encode = restrict(pos, start_rec, sleep_times);
pos_recall = restrict(pos, end_sleep_times, pos.tvec(end));
pos_sleep = restrict(pos, sleep_times, end_sleep_times);

pos_encode_trials = restrict(pos, t_start_encode+60, t_end_encode);
pos_encode_iti = restrict(pos, t_start_encode, t_start_encode+60);

pos_recall_trials = restrict(pos, t_start_recall+60, t_end_recall);
pos_recall_iti = restrict(pos, t_start_recall, t_start_recall+60);

fprintf('Encoding dur: %.0fmins \nSleep dur: %.0fmins \nRecall dur: %.0fmins\n', (pos_encode.tvec(end) - pos_encode.tvec(1))/60,  (pos_sleep.tvec(end) - pos_sleep.tvec(1))/60, (pos_recall.tvec(end) - pos_recall.tvec(1))/60)

clear pos

% test plot
figure(101)
clf
subplot(3,3,1)
hold on
plot(pos_encode_trials.data(1,:), pos_encode_trials.data(2,:), '.', 'color', c_ord(1,:));
plot(pos_encode_iti.data(1,:), pos_encode_iti.data(2,:), '.', 'color', c_ord(2,:));
legend({'Trials'; 'ITI'})
title('Ecoding phase')
subplot(3,3,2)
plot(pos_sleep.data(1,:), pos_sleep.data(2,:), '.k');
title('Sleep')

subplot(3,3,3)
hold on
plot(pos_recall_trials.data(1,:), pos_recall_trials.data(2,:), '.', 'color', c_ord(1,:));
plot(pos_recall_iti.data(1,:), pos_recall_iti.data(2,:), '.', 'color', c_ord(2,:));
title('Recall phase')

%% get some tvec info
cd(raw_data_dir)
tvec = readNPY('timestamps.npy');


%% load some spikes and convert them to the TS format.
cd(spike_dir)

S = OE_phy2TS;

% convert S.t to time using timestamps.
for ii = 1:length(S.t)
    S.t{ii} = tvec(S.t{ii});
    
end

%restrict to encoding/sleep/recall
S_encode = restrict(S, start_rec, sleep_times);
S_recall = restrict(S, end_sleep_times, end_rec);
S_sleep = restrict(S, sleep_times, end_sleep_times);

S_encode_trials = restrict(S, t_start_encode+60, t_end_encode);
S_encode_iti = restrict(S, t_start_encode, t_start_encode+60);

S_recall_trials = restrict(S, t_start_recall+60, t_end_recall);
S_recall_iti = restrict(S, t_start_recall, t_start_recall+60);

%% plot the spike rasters

figure(101)
subplot(3,3,7:9)
plot(S);
% for ii = 1:length(REM_sleep_times)
% xline([REM_sleep_times(ii)], 'color', c_ord(2,:))
% end
xline([sleep_times], 'color', c_ord(3,:), 'linewidth', 5)
xline([end_sleep_times], 'color', c_ord(3,:), 'linewidth', 5)



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
SET_xmin = 0; SET_ymin = 0; % set up bins
SET_xmax = 500; SET_ymax = 700;
SET_xBinSz = 20; SET_yBinSz = 20;

 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([1 1],1); % 2d gaussian in bins

% compute occupancy encode
occ_hist_encode = histcn(pos_encode_trials.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist_encode = conv2(occ_hist_encode,kernel,'same');
 
no_occ_idx_encode = find(occ_hist_encode < 7.5); % NaN out bins never visited
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
        figure(102+ii)
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

cell_id = 549; 
ii = 48; 
ii = 33

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
    
    