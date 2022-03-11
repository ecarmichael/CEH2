%% sandbox_convert MAze to 1D

addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

cd('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3');

load('behav_DLC.mat');
load('ms_trk.mat'); 
load('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\Common_CoorD.mat')

%% office computer

addpath(genpath('/home/williamslab/Documents/Github/CEH2'));
addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));

cd('/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/Maze_Ca/inter/pv1254/2021_12_17_pv1254_MZD3');

load('behav_DLC.mat');
load('ms_trk.mat'); 



%% remove inactive cells

keep_idx = sum(ms_trk.Binary, 1) >0; 

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv'; 
ms_trk_act = MS_Remove_trace(cfg_rem, ms_trk); 


%% play with linerization using one direction at a time


pos = MS_behav2tsd(behav); 
pos.tvec = pos.tvec/1000; 
% limit to movement
linspeed = getLinSpd([],pos); % linear speed

linspeed.data = smooth(linspeed.data, floor(1/mode(diff(linspeed.tvec)))*.5);
move_idx = linspeed.data > 5; 


pos_move  = pos;
pos_move.tvec = pos_move.tvec(move_idx);
pos_move.data = pos_move.data(:, move_idx); 

ms_trk_move = ms_trk_act; 
ms_trk_move.time = ms_trk_act.time(move_idx);
ms_trk_move.Binary = ms_trk_act.Binary(move_idx,:);

% convert ms_trk_move to tsd
ms_tsd = tsd(ms_trk_move.time'/1000, ms_trk.Binary', 'Binary'); 

%% split L and R trials
trials_R_s = sort([behav.trials_ms.start_CL; behav.trials_ms.start_FL])/1000;
trials_R_e = sort([behav.trials_ms.rew_L])/1000;

trials_L_s = sort([behav.trials_ms.start_CR; behav.trials_ms.start_FR])/1000;
trials_L_e = sort([behav.trials_ms.rew_R])/1000;


pos_L = restrict(pos_move, trials_L_s, trials_L_e); 
pos_R = restrict(pos_move, trials_R_s, trials_R_e); 

ms_L = restrict(ms_tsd, trials_L_s, trials_L_e);

ms_R = restrict(ms_tsd, trials_R_s, trials_R_e);

figure(101)
hold on
plot(pos_L.data(1,:), pos_L.data(2,:), 'b');
plot(pos_R.data(1,:), pos_R.data(2,:), 'r');
legend({'left', 'right'})

%% exclude times when they are not on the track or errors
L_limits = [55 24] ; % x and y limits
R_limits = [45 24] ; % x and y limits

keep_idx = (pos_L.data(1,:) <= L_limits(1)) & (pos_L.data(2,:) >= L_limits(2));

pos_L_cut = pos_L;
pos_L_cut.data = pos_L.data(:, keep_idx); 
pos_L_cut.tvec = pos_L.tvec(keep_idx); 

ms_L_cut = ms_L;
ms_L_cut.data = ms_L.data(:, keep_idx);
ms_L_cut.tvec = ms_L.tvec(:,keep_idx); 


keep_idx = (pos_R.data(1,:) >= R_limits(1)) & (pos_R.data(2,:) >= R_limits(2));

pos_R_cut = pos_R;
pos_R_cut.data = pos_R.data(:, keep_idx); 
pos_R_cut.tvec = pos_R.tvec(keep_idx); 

ms_R_cut = ms_R;
ms_R_cut.data = ms_R.data(:, keep_idx);
ms_R_cut.tvec = ms_R.tvec(:,keep_idx); 

figure(101)
clf
subplot(2,2,1)
hold on
plot(pos_L.data(1,:), pos_L.data(2,:), '.r')
plot(pos_R.data(1,:), pos_R.data(2,:), '.b')
xlim([min([pos_L.data(1,:), pos_move.data(1,:)]), max([pos_L.data(1,:), pos_move.data(1,:)])]); 
ylim([min([pos_L.data(2,:), pos_move.data(2,:)]), max([pos_L.data(2,:), pos_move.data(2,:)])]); 

legend({'left', 'right'})

subplot(2,2,2)
hold on
plot(pos_L_cut.data(1,:), pos_L_cut.data(2,:), '.r')
plot(pos_R_cut.data(1,:), pos_R_cut.data(2,:), '.b')
xlim([min([pos_L.data(1,:), pos_move.data(1,:)]), max([pos_L.data(1,:), pos_move.data(1,:)])]); 
ylim([min([pos_L.data(2,:), pos_move.data(2,:)]), max([pos_L.data(2,:), pos_move.data(2,:)])]); 
legend({'left', 'right'})

%% snap position to ideal

P = behav.CoorD_L.coord'; 
PQ = pos_L_cut.data'; 
[k, ~] = dsearchn(P, PQ);

pos_snap_L = pos_L_cut; 
pos_snap_L.data(1,:) = P(k, 1);
pos_snap_L.data(2,:) = P(k, 2);

P = behav.CoorD_R.coord'; 
PQ = pos_R_cut.data'; 
[k, ~] = dsearchn(P, PQ);

pos_snap_R = pos_R_cut; 
pos_snap_R.data(1,:) = P(k, 1);
pos_snap_R.data(2,:) = P(k, 2);



figure(101)
subplot(2,2,2)
% hold on
% plot(P(:,1),P(:,2),'ko')
% hold on
% plot(PQ(:,1),PQ(:,2),'*g')
hold on
plot(pos_snap_L.data(1,:),pos_snap_L.data(2,:),'ok', 'markersize', 10, 'markeredgecolor', 'k')
plot(pos_snap_L.data(1,:),pos_snap_L.data(2,:),'xm', 'markersize', 10)
xlim([min([pos_L.data(1,:), pos_move.data(1,:)]), max([pos_L.data(1,:), pos_move.data(1,:)])]); 
ylim([min([pos_L.data(2,:), pos_move.data(2,:)]), max([pos_L.data(2,:), pos_move.data(2,:)])]); 
legend('Nearest Pos L')



% hold on
% plot(P(:,1),P(:,2),'ko')
% hold on
% plot(PQ(:,1),PQ(:,2),'*g')
hold on
plot(pos_snap_R.data(1,:),pos_snap_R.data(2,:),'ok', 'markersize', 10, 'markeredgecolor', 'k')
plot(pos_snap_R.data(1,:),pos_snap_R.data(2,:),'xc', 'markersize', 10)
xlim([min([pos_L.data(1,:), pos_move.data(1,:)]), max([pos_L.data(1,:), pos_move.data(1,:)])]); 
ylim([min([pos_L.data(2,:), pos_move.data(2,:)]), max([pos_L.data(2,:), pos_move.data(2,:)])]); 
legend('Nearest Pos R')

%% standardize and linearize
% load(
cfg = [];
cfg.binsize = .5;
cfg.run_dist = Common_CoorD.run_dist;
SCoorD_R = StandardizeCoord(cfg, behav.CoorD_L, Common_CoorD.run_dist);
SCoorD_L = StandardizeCoord(cfg, behav.CoorD_R, Common_CoorD.run_dist);

%
cfg = [];
cfg.Coord = SCoorD_L.coord;
linpos.L = LinearizePos(cfg,pos_L_cut, SCoorD_L);
cfg.Coord = SCoorD_R.coord;
linpos.R = LinearizePos(cfg,pos_R_cut, SCoorD_R);

figure(1)
subplot(2,2,3)
plot(linpos.L.tvec, linpos.L.data, '.r')
subplot(2,2,4)
plot(linpos.R.tvec, linpos.R.data, '.b')

%% compute tuning curves for each cell
bin_s = 1; 
bins = floor(min(linpos.L.data)):bin_s:ceil(max(linpos.L.data));

xbins = floor(min(pos_L_cut.data(1,:))):bin_s:ceil(max(pos_L_cut.data(1,:)))
ybins = floor(min(pos_L_cut.data(2,:))):bin_s:ceil(max(pos_L_cut.data(2,:)))
nShuff = 500;

MI = []; posterior = []; occupancy = []; p_active = []; 
for iCell = size(ms_L_cut.data,1):-1:1
    
    if sum(ms_L_cut.data(iCell, :)) ==0
        MI(iCell) = 0;
        posterior(iCell,:) = NaN(iCell, length(bins));
        TC(iCell, :) = posterior(iCell,:);
        p_active(iCell) = 0;
        
    else
%         [MI(iCell), posterior(iCell,:),occupancy,p_active(iCell), TC(iCell,:)] = MS_get_spatial_information(ms_L_cut.data(iCell, :), linpos.L.data, bins);
        [MI(iCell), posterior(iCell,:),occupancy,p_active(iCell), TC(iCell,:)] = GE_extract_1D_information(logical(ms_L_cut.data(iCell, :))', linpos.L.data', bins', ones(1, length(ms_L_cut.data(iCell, :))));

        % convert likel
        for ishuff = nShuff:-1:1
            shuf_data = circshift(ms_L_cut.data(iCell, :), round(MS_randn_range(1,1,1,length(ms_L_cut.data(iCell, :)))));
            
            [shuff_MI(ishuff), ~, ~, ~, shuff_TC(ishuff,:)] = MS_get_spatial_information(shuf_data, linpos.L.data, bins);
            
            pval = sum(shuff_TC(ishuff,:) > TC(iCell,:),1)/nShuff;
            
            sig_TC(iCell,:) = TC(iCell,:); 
            
            sig_TC(iCell,pval > 0.05) = 0; 
            
%             if sum(sig_TC < 0.05)
            sig_idx = 1; 
            
            % test plot if needed
            
            figure(102)
            imagesc(xbins, ybins, o')
            hold on
            plot(pos_snap_L.data(1,:), pos_snap_L.data(2,:), '.');
            set(gca, 'YDir', 'reverse')
            
            
        end
        
    end
end


end