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

%% play with linerization using one direction at a time

% limit to movement
linspeed = getLinSpd([],pos); % linear speed

linspeed.data = smooth(linspeed.data, floor(1/mode(diff(linspeed.tvec)))*.5);
move_idx = linspeed.data > 5; 


pos_move  = pos;
pos_move.tvec = pos_move.tvec(move_idx);
pos_move.data = pos_move.data(:, move_idx); 

spike_x = interp1(pos.tvec, pos.data(1,:), S.t{1}, 'linear'); 
spike_y = interp1(pos.tvec, pos.data(2,:), S.t{1}, 'linear'); 

figure(101)
subplot(2,2,1)
hold on
plot(pos.data(1,:), pos.data(2,:), '.', 'color', [.2 .2 .2]); 
plot(spike_x, spike_y, '.r')

t_types = unique(maze.trials.types); 
c_trial = linspecer(length(t_types)); 
for ii = 1:length(maze.trials.types)
    this_pos = restrict(pos, maze.trials.times(ii,1), maze.trials.times(ii, 2)); 
    idx = find(contains(t_types, maze.trials.types{ii}));
    plot(this_pos.data(1,:), this_pos.data(2,:), '.', 'color', c_trial(idx,:))
    
end

%% split L and R trials
R_idx = contains( maze.trials.types, {'FR correct'; 'CR correct'; 'CL error'}); 
L_idx = contains( maze.trials.types, {'FL correct'; 'CL correct'; 'CR error'}); 
R_iv = iv(maze.trials.times(R_idx,1), maze.trials.times(R_idx,2));
L_iv = iv(maze.trials.times(L_idx,1), maze.trials.times(L_idx,2));

pos_L = restrict(pos, L_iv); 
pos_R = restrict(pos, R_iv); 

S_L = restrict(S, L_iv);

S_R = restrict(S, R_iv);

figure(101)
subplot(2,2,2)
hold on
plot(pos_L.data(1,:), pos_L.data(2,:), '.b');
plot(pos_R.data(1,:), pos_R.data(2,:), '.r');
legend({'left', 'right'})


%% snap position to ideal

PL = Common_CoorD.CoorD_L.coord'; 
PQL = pos_L.data(1:2,:)'; 
[k, ~] = dsearchn(PL, PQL);

pos_snap_L = pos_L; 
pos_snap_L.data(1,:) = PL(k, 1);
pos_snap_L.data(2,:) = PL(k, 2);

PR = Common_CoorD.CoorD_R.coord'; 
PQR = pos_R.data(1:2,:)'; 
[k, ~] = dsearchn(PR, PQR);

pos_snap_R = pos_R; 
pos_snap_R.data(1,:) = PR(k, 1);
pos_snap_R.data(2,:) = PR(k, 2);



figure(101)
subplot(2,2,2)
cla
hold on
% plot(PL(:,1),PL(:,2),'ko')
% hold on
plot(PQL(:,1),PQL(:,2),'*b', 'markersize', 2)

% hold on
hl1 = plot(pos_snap_L.data(1,:),pos_snap_L.data(2,:),'ok', 'markersize', 10, 'markeredgecolor', 'k');
hl2 = plot(pos_snap_L.data(1,:),pos_snap_L.data(2,:),'xc', 'markersize', 10);
xlim([min([pos_L.data(1,:), pos_move.data(1,:)]), max([pos_L.data(1,:), pos_move.data(1,:)])]); 
ylim([min([pos_L.data(2,:), pos_move.data(2,:)]), max([pos_L.data(2,:), pos_move.data(2,:)])]); 

x_lim = xlim; 
y_lim = ylim; 

hold on
plot(PQR(:,1),PQR(:,2),'*r', 'markersize', 2)
hr1 = plot(pos_snap_R.data(1,:),pos_snap_R.data(2,:),'dk', 'markersize', 10, 'markeredgecolor', 'k');
hr2 = plot(pos_snap_R.data(1,:),pos_snap_R.data(2,:),'xm', 'markersize', 10);
xlim([min([pos_R.data(1,:), pos_move.data(1,:)]), max([pos_R.data(1,:), pos_move.data(1,:)])]); 
ylim([min([pos_R.data(2,:), pos_move.data(2,:)]), max([pos_R.data(2,:), pos_move.data(2,:)])]); 
legend([hl1 hr1, hl2, hr2],{'left', 'right','Nearest Pos L', 'Nearest Pos R'}, 'location', 'southeast')


%% standardize and linearize
% load(
cfg = [];
cfg.binsize = .5;
cfg.run_dist = Common_CoorD.run_dist;
SCoorD_R = StandardizeCoord(cfg, Common_CoorD.CoorD_L, Common_CoorD.run_dist);
SCoorD_L = StandardizeCoord(cfg, Common_CoorD.CoorD_R, Common_CoorD.run_dist);

%
cfg = [];
cfg.Coord = SCoorD_L.coord;
linpos.L = LinearizePos(cfg,pos_L, SCoorD_L);
cfg.Coord = SCoorD_R.coord;
linpos.R = LinearizePos(cfg,pos_R, SCoorD_R);

subplot(2,2,3)
cla
hold on 
plot(SCoorD_L.coord(1,:), SCoorD_L.coord(2,:), '.b');
plot(SCoorD_R.coord(1,:), SCoorD_R.coord(2,:), '.r');
xlim(x_lim);
ylim(y_lim);


subplot(2,2,4)
hold on
hlinl = plot(linpos.L.tvec, linpos.L.data, '.b');
xlim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])])
hlinr = plot(linpos.R.tvec, linpos.R.data, '.r');
xlim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])]);
legend([hlinl, hlinr], 'Linearized left trajectories', 'Linearized right trajectories'); 

%% compute tuning curves for each cell
bin_s = 3; 
bins = floor(min([linpos.R.data, linpos.L.data])):bin_s:39%ceil(max([linpos.R.data, linpos.L.data]));

figure(102)
subplot(2,2,1)
imagesc(bins,1:2, [histc(linpos.L.data, bins)/mode(diff(linpos.L.tvec)); histc(linpos.R.data, bins)/mode(diff(linpos.R.tvec))]/60); 
set(gca,'Ytick', 1:2, 'YTickLabel', {'left occupancy', 'right occupancy'}); 
c = colorbar; 
% c.Label = 'time (s)'; 

% xbins = floor(min(pos_L_cut.data(1,:))):bin_s:ceil(max(pos_L_cut.data(1,:)))
% ybins = floor(min(pos_L_cut.data(2,:))):bin_s:ceil(max(pos_L_cut.data(2,:)))
nShuff = 500;

MI = []; posterior = []; occupancy = []; p_active = []; TC = []; 
for iCell = size(ms_L_cut.data,1):-1:1
    
    if sum(ms_L_cut.data(iCell, :)) ==0
        MI(iCell) = 0;
        posterior(iCell,:) = NaN(1, length(bins)-1);
        TC(iCell, :) = posterior(iCell,:);
        p_active(iCell) = 0;
        
    else
%         [MI(iCell), posterior(iCell,:),occupancy,p_active(iCell), TC(iCell,:)] = MS_get_spatial_information(ms_L_cut.data(iCell, :), linpos.L.data, bins);
        [MI(iCell), posterior(iCell,:),occupancy,p_active(iCell), TC(iCell,:)] = GE_extract_1D_information((ms_R_cut.data(iCell, :)>0)', linpos.R.data', bins', ones(length(ms_R_cut.data(iCell, :)),1));

        figure(191)
        hold on
        
        
        
        % convert likel
        for ishuff = nShuff:-1:1
            shuf_data = circshift(ms_L_cut.data(iCell, :), round(MS_randn_range(1,1,1,length(ms_L_cut.data(iCell, :)))));
            
            [shuff_MI(ishuff), ~, ~, ~, shuff_TC(ishuff,:)] = GE_extract_1D_information((shuf_data>0)', linpos.L.data', bins',ones(length(ms_L_cut.data(iCell, :)),1));
            
              
            % test plot if needed
        end
        
        figure(1010)
        clf
        hold on
        plot(bins(1:end-1) + diff(bins), TC(iCell,:), 'b');
        plot(bins(1:end-1) + diff(bins), mean(shuff_TC,1), '--r')
        
        
        pval = sum(shuff_TC(ishuff,:) > TC(iCell,:),1)/nShuff;
        
        sig_TC(iCell,:) = TC(iCell,:);
        
        sig_TC(iCell,pval > 0.05) = 0;
        
        if sum(sig_TC < 0.05)
            sig_idx = 1;
    end
        
            figure(102)
            imagesc(xbins, ybins, occupancy')
            hold on
            plot(pos_snap_L.data(1,:), pos_snap_L.data(2,:), '.');
            set(gca, 'YDir', 'reverse')

    end
end


