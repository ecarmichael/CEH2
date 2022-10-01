%% dsub_maze_1d

addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

cd('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\dSubiculum\inProcess\MD3\MD3_2022_07_26_D27');

%% load stuff
evts = LoadEvents([]); 

blocks = MS_get_waze_blocks(evts); 

cfg_pos = [];
cfg_pos.convFact = [6.95 6.95];
pos = MS_LoadPos(cfg_pos); 
pos.data  = pos.data(1:2,:); 
pos.label = pos.label(1:2); 
load('maze.mat'); 

meta = MS_Load_meta();
load('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\dSubiculum\inProcess\MD3\Common_CoorD.mat')

cfg = [];
cfg.fc = {'TT2_01_3.t'}; 
cfg.getTTnumbers = 0; 
S = LoadSpikes(cfg); 


pos = restrict(pos, blocks.W_maze(1), blocks.W_maze(2)); 
S = restrict(S, blocks.W_maze(1), blocks.W_maze(2)); 
evts = restrict(evts, blocks.W_maze(1), blocks.W_maze(2)); 

%%
linspeed = getLinSpd([],pos); % linear speed
 
% Threshold speed
cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 1; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg,linspeed); % only keep intervals with speed above thresh
 
% Restrict data so it includes fast intervals only
pos_fast = restrict(pos,iv_fast);
S_fast = restrict(S,iv_fast);

%%


spike_x = interp1(pos_fast.tvec, pos_fast.data(1,:), S_fast.t{1}, 'linear'); 
spike_y = interp1(pos_fast.tvec, pos_fast.data(2,:), S_fast.t{1}, 'linear'); 

figure(101)
subplot(2,2,1)
hold on
plot(pos_fast.data(1,:), pos_fast.data(2,:), '.', 'color', [.2 .2 .2]); 
plot(spike_x, spike_y, '.r')
%% split L and R trials
R_idx = contains( maze.trials.types, {'FR correct'; 'CR correct'; 'CL error'}); 
L_idx = contains( maze.trials.types, {'FL correct'; 'CL correct'; 'CR error'}); 
R_iv = iv(maze.trials.times(R_idx,1), maze.trials.times(R_idx,2));
L_iv = iv(maze.trials.times(L_idx,1), maze.trials.times(L_idx,2));

pos_L = restrict(pos_fast, L_iv); 
pos_R = restrict(pos_fast, R_iv); 

S_L = restrict(S_fast, L_iv);
S_R = restrict(S_fast, R_iv);

figure(101)
subplot(2,2,2)
hold on
plot(pos_L.data(1,:), pos_L.data(2,:), '.b');
plot(pos_R.data(1,:), pos_R.data(2,:), '.r');
legend({'left', 'right'})
%% standardize and linearize
cfg = [];
cfg.binsize = .5;
cfg.run_dist = Common_CoorD.run_dist;
coordL_cm = []; coordL_cm.coord = Common_CoorD.CoorD_L.coord; coordL_cm.units = 'cm'; 
coordR_cm = []; coordR_cm.coord = Common_CoorD.CoorD_R.coord;  coordR_cm.units = 'cm'; 

coord.L = StandardizeCoord(cfg,coordL_cm, Common_CoorD.run_dist);
coord.R = StandardizeCoord(cfg,coordR_cm, Common_CoorD.run_dist);


%
cfg = [];
cfg.Coord = coord.L;
linpos.L = LinearizePos([],pos_L, coord.L);
cfg.Coord = coord.R;
linpos.R = LinearizePos([],pos_R, coord.R);

subplot(2,2,3)
cla
hold on 
hl = plot(coordL_cm.coord(1,:), coordL_cm.coord(2,:), '.b');
h2 = plot(coordR_cm.coord(1,:), coordR_cm.coord(2,:), '.r');

% xlim(x_lim);
% ylim(y_lim);

spk_idx_L = nearest_idx3(S_L.t{1}, linpos.L.tvec); 
spk_idx_R = nearest_idx3(S_R.t{1}, linpos.R.tvec); 

subplot(2,2,4)
c_ord = linspecer(6); 
cla
hold on
rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 0,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 25], 'facecolor', [c_ord(2,:), .2], 'edgecolor', 'none')
rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 25,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 5], 'facecolor', [c_ord(6,:), .2], 'edgecolor', 'none')

rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 30,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 20], 'facecolor', [c_ord(3,:), .2], 'edgecolor', 'none')
rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 50,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 20], 'facecolor', [c_ord(5,:), .2], 'edgecolor', 'none')
rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 70,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 20], 'facecolor', [c_ord(4,:), .2], 'edgecolor', 'none')
rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 90,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 10], 'facecolor', [c_ord(1,:), .2], 'edgecolor', 'none')

hlinl = plot(linpos.L.tvec, linpos.L.data, '.', 'color', [0.2 0.2 .2]);
hlinls = plot(linpos.L.tvec(spk_idx_L), linpos.L.data(spk_idx_L), 'x', 'color', [0 0 1]);

xlim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])])
hlinr = plot(linpos.R.tvec, linpos.R.data, '.', 'color', [.2 0.2 0.2]);
hlinrs = plot(linpos.R.tvec(spk_idx_R), linpos.R.data(spk_idx_R), 'x', 'color', [1 0 0]);

xlim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])]);
legend([hlinl, hlinr, hlinls, hlinrs], 'Linearized left trajectories', 'Linearized right trajectories', 'Left spikes', 'Right spikes'); 
yyaxis right
cla
ylim([0 100])
set(gca, 'ytick', [12.5 27.5 40 60 80 95], 'YTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})

%% generate TCs

% make spatial bins
bin_s = 2.5; 
bins = 0-(bin_s/2):bin_s:100+(bin_s/2);
bin_c = bins(1:end-1)+bin_s/2;

% gk cfgs
gk_win = 5;
gk_sd = 2.5; 
gkernal = gausskernel(gk_win/bin_s, gk_sd/bin_s); 

Fs = 30; 
occ_L = histc(linpos.L.data, bins)/Fs; 
occ_L = trim_histc(occ_L); 
occ_L(occ_L==0) = NaN;
occ_L = conv(occ_L, gkernal, 'same'); 

occ_R = histc(linpos.R.data, bins)/Fs; 
occ_R = trim_histc(occ_R); 
occ_R(occ_R==0) = NaN;
occ_R = conv(occ_R, gkernal, 'same'); 

% spk_hist
spk_L = interp1(linpos.L.tvec, linpos.L.data, S_L.t{1}, 'linear');
spk_h_L = histc(spk_L', bins); 
spk_h_L = conv(spk_h_L(1:end-1), gkernal, 'same');  
spk_tc_L = spk_h_L./occ_L; 

spk_R = interp1(linpos.R.tvec, linpos.R.data, S_R.t{1}, 'linear');
spk_h_R = histc(spk_R', bins); 
spk_h_R = conv(spk_h_R(1:end-1), gkernal, 'same');  
spk_tc_R = spk_h_R./occ_R; 


%% plot TCs
figure(102)
subplot(6,1,1:4)
c_ord = linspecer(6); 
cla
hold on
tvec_L = linpos.L.tvec - linpos.L.tvec(1); 
tvec_R = linpos.R.tvec - linpos.R.tvec(1); 

rectangle('position', [0, min([tvec_L, tvec_R]), 25,  max([tvec_L, tvec_R]) - min([tvec_L, tvec_R])], 'facecolor', [c_ord(2,:), .2], 'edgecolor', 'none')
rectangle('position', [25, min([tvec_L, tvec_R]), 5,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec])], 'facecolor', [c_ord(6,:), .2], 'edgecolor', 'none')

rectangle('position', [30, min([tvec_L, tvec_R]), 20, max([tvec_L, tvec_R]) - min([tvec_L, tvec_R])], 'facecolor', [c_ord(3,:), .2], 'edgecolor', 'none')
rectangle('position', [50, min([tvec_L, tvec_R]), 20, max([tvec_L, tvec_R]) - min([tvec_L, tvec_R])], 'facecolor', [c_ord(5,:), .2], 'edgecolor', 'none')
rectangle('position', [70, min([tvec_L, tvec_R]), 20, max([tvec_L, tvec_R]) - min([tvec_L, tvec_R])], 'facecolor', [c_ord(4,:), .2], 'edgecolor', 'none')
rectangle('position', [90, min([tvec_L, tvec_R]), 10, max([tvec_L, tvec_R]) - min([tvec_L, tvec_R])], 'facecolor', [c_ord(1,:), .2], 'edgecolor', 'none')


hlinl = plot(linpos.L.data, tvec_L, '.', 'color', [0.7 0.7 .7]);
hlinls = plot(linpos.L.data(spk_idx_L), tvec_L(spk_idx_L), 'x', 'color', [0 0 1]);

hlinr = plot(linpos.R.data, tvec_R, '.', 'color', [.7 0.7 0.7]);
hlinrs = plot(linpos.R.data(spk_idx_R),tvec_R(spk_idx_R),  'x', 'color', [1 0 0]);
ylim([min([tvec_L, tvec_R]), max([tvec_L, tvec_R])])

% ylim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])]);
set(gca, 'YDir', 'reverse')
ylabel('time (s)')
legend([hlinl, hlinr, hlinls, hlinrs], 'Linearized left trajectories', 'Linearized right trajectories', 'Left spikes', 'Right spikes'); 

set(gca, 'xtick', [12.5 27.5 40 60 80 95], 'xTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})

subplot(6,1,5)
cla
hold on
yyaxis left
plot(bin_c, occ_L,  '-', 'color', [.7 .7 .7])
ylim([min([occ_L occ_R]), max([occ_L occ_R])]); 
set(gca, 'YColor', 'k')
ylabel('Occupancy L (s)')
yyaxis right
plot(bin_c,spk_tc_L, '-', 'color', c_ord(2,:))
ylim([min([spk_tc_L spk_tc_R]), max([spk_tc_L spk_tc_R])]); 
ylabel('firing rate (Hz)')
set(gca, 'YColor', c_ord(1,:))


subplot(6,1,6)
cla
hold on
yyaxis left
plot(bin_c, occ_R,  '-', 'color', [.7 .7 .7])
ylim([min([occ_L occ_R]), max([occ_L occ_R])]); 
set(gca, 'YColor', 'k')
ylabel('Occupancy R (s)')
yyaxis right
plot(bin_c,spk_tc_R, '-', 'color', c_ord(1,:))
ylim([min([spk_tc_L spk_tc_R]), max([spk_tc_L spk_tc_R])]); 
ylabel('firing rate (Hz)')
set(gca, 'YColor', c_ord(1,:))

set(gca, 'xtick', [12.5 27.5 40 60 80 95], 'xTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})

% get a rate map for each lap. 
