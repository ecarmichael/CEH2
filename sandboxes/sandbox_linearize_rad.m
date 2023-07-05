%% sandbox Linearize radial maze

cd('/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Radial_ephys/M29/M29_2023-04-14_12-13-30_Rad10/M29_2023-04-14_Rad10')

pos = LoadPos([]);

evt = LoadEvents([]);


start_rec = find(contains(evt.label, 'Starting Recording'));
stop_rec = find(contains(evt.label, 'Stopping Recording'));
sleep = find(contains(evt.label, 'sleep'));




pos_encode = restrict(pos, evt.t{start_rec},  evt.t{sleep}(1));
pos_recall= restrict(pos,  evt.t{sleep}(2) , evt.t{stop_rec});


%%
% set up bins

SET_xmin = 0; SET_ymin = 0; % set up bins
SET_xmax = 500; SET_ymax = 700;
SET_xBinSz = 20; SET_yBinSz = 20;

 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([1 1],1); % 2d gaussian in bins

% compute occupancy encode
occ_hist_encode = histcn(pos_encode.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist_encode = conv2(occ_hist_encode,kernel,'same');
 
no_occ_idx_encode = find(occ_hist_encode < 7.5); % NaN out bins never visited
occ_hist_encode(no_occ_idx_encode) = NaN;
 
occ_hist_encode = occ_hist_encode .* (1/30); % convert samples to seconds using video frame rate (30 Hz)

% compute occupancy recall
occ_hist_recall = histcn(pos_recall.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist_recall = conv2(occ_hist_recall,kernel,'same');

no_occ_idx_recall = find(occ_hist_recall < 7.5); % NaN out bins visited for less than 1 s
occ_hist_recall(no_occ_idx_recall) = NaN;
 
occ_hist_recall = occ_hist_recall .* (1/30); % convert samples to seconds using video frame rate (30 Hz)



%% plot check


figure(110)
clf
subplot(2,2,1)
plot(pos_encode.data(1,:), pos_encode.data(2,:), '.b')


subplot(2,2,2)
plot(pos_recall.data(1,:), pos_recall.data(2,:), '.r')

subplot(2,2,3)
pcolor(occ_hist_encode)


subplot(2,2,4)
pcolor(occ_hist_recall)
