function TC_out = dSub_gen_1d_TC(cfg_in, data_dir, Common_CoorD)
%% MS_gen_1d_TC: generate a 1d tuning curve by linea
%
%
%
%    Inputs:
%    -
%
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2022-10-01   initial version
%
%
%
%% initialize
if nargin < 1
    cfg_in = [];
    data_dir = cd;
    Common_CoorD = [];
elseif nargin ==2
    data_dir = cd;
    Common_CoorD = [];
elseif nargin < 3
    Common_CoorD = [];
end

cd(data_dir)

% configs
cfg_def = [];
cfg_def.bin_s = 2.5;
cfg_def.gk_win = 5;
cfg_def.gk_sd = 2.5;
cfg_def.plot = 1; % plot flag.

cfg = ProcessConfig(cfg_def, cfg_in);


% make spatial bins
bins = 0-(cfg.bin_s/2):cfg.bin_s:100+(cfg.bin_s/2);
bin_c = bins(1:end-1)+cfg.bin_s/2;

% gk cfgs
gkernal = gausskernel(cfg.gk_win/cfg.bin_s, cfg.gk_sd/cfg.bin_s);

%% load stuff

fname = strsplit(data_dir, filesep);
fname = strrep(fname{end}, '-', '_');
parts = strsplit(fname, '_'); 
sess = parts{end}; 

evts = LoadEvents([]);

blocks = MS_get_waze_blocks(evts);

cfg_pos = [];
cfg_pos.convFact = [6.95 6.95];
pos = MS_LoadPos(cfg_pos);
pos.data  = pos.data(1:2,:);
pos.label = pos.label(1:2);
load('maze.mat');

% meta = MS_Load_meta();
if isempty(Common_CoorD)
    load('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\dSubiculum\inProcess\MD3\Common_CoorD.mat')
end

if str2double(sess(2:end)) < 14
    Common_CoorD.CoorD_L.coord(1,:) = Common_CoorD.CoorD_L.coord(1,:) - 4; 
    Common_CoorD.CoorD_R.coord(1,:) = Common_CoorD.CoorD_R.coord(1,:) - 4; 
end

cfg_S = [];
% cfg.fc = {'TT2_01_3.t'};
cfg_S.getTTnumbers = 0;
% cfg.min_cluster_quality = 4;
S = LoadSpikes(cfg_S);

% remove low quality cells
for iS = length(S.t):-1:1
    if strcmp(S.label{iS}(end-2), '5')
        fprintf('Removing cell: %s due to low quality (%s)\n', S.label{iS}, S.label{iS}(end-2))
        S.t(iS) = [];
        S.label(iS) = [];
    end
end



pos = restrict(pos, blocks.W_maze(1), blocks.W_maze(2));
S = restrict(S, blocks.W_maze(1), blocks.W_maze(2));
evts = restrict(evts, blocks.W_maze(1), blocks.W_maze(2));

%%
linspeed = getLinSpd([],pos); % linear speed

% Threshold speed
cfg_l = []; cfg_l.method = 'raw'; cfg_l.operation = '>'; cfg_l.threshold = 1; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg_l,linspeed); % only keep intervals with speed above thresh

% Restrict data so it includes fast intervals only
pos_fast = restrict(pos,iv_fast);
S_fast = restrict(S,iv_fast);

%%
%
%
% spike_x = interp1(pos_fast.tvec, pos_fast.data(1,:), S_fast.t{1}, 'linear');
% spike_y = interp1(pos_fast.tvec, pos_fast.data(2,:), S_fast.t{1}, 'linear');
%
% figure(101)
% subplot(2,2,1)
% hold on
% plot(pos_fast.data(1,:), pos_fast.data(2,:), '.', 'color', [.2 .2 .2]);
% plot(spike_x, spike_y, '.r')
%% split L and R trials
R_idx = contains( maze.trials.types, {'FR correct'; 'CR correct'; 'CL error'});
L_idx = contains( maze.trials.types, {'FL correct'; 'CL correct'; 'CR error'});
R_iv = iv(maze.trials.times(R_idx,1), maze.trials.times(R_idx,2)+5);
L_iv = iv(maze.trials.times(L_idx,1), maze.trials.times(L_idx,2)+5);

pos_L = restrict(pos_fast, L_iv);
pos_R = restrict(pos_fast, R_iv);

S_L = restrict(S_fast, L_iv);
S_R = restrict(S_fast, R_iv);

% figure(101)
% clf
% subplot(2,2,2)
% hold on
% plot(pos_L.data(1,:), pos_L.data(2,:), '.b');
% plot(pos_R.data(1,:), pos_R.data(2,:), '.r');
% legend({'left', 'right'})
%% standardize and linearize
cfg_c = [];
cfg_c.binsize = .5;
cfg_c.run_dist = Common_CoorD.run_dist;
coordL_cm = []; coordL_cm.coord = Common_CoorD.CoorD_L.coord; coordL_cm.units = 'cm';
coordR_cm = []; coordR_cm.coord = Common_CoorD.CoorD_R.coord;  coordR_cm.units = 'cm';

coord.L = StandardizeCoord(cfg_c,coordL_cm, Common_CoorD.run_dist);
coord.R = StandardizeCoord(cfg_c,coordR_cm, Common_CoorD.run_dist);


% %
% cfg_l = [];
% cfg_l.Coord = coord.L;
linpos.L = LinearizePos([],pos_L, coord.L);
% cfg_l.Coord = coord.R;
linpos.R = LinearizePos([],pos_R, coord.R);

% subplot(2,2,3)
% cla
% hold on
% hl = plot(coordL_cm.coord(1,:), coordL_cm.coord(2,:), '.b');
% h2 = plot(coordR_cm.coord(1,:), coordR_cm.coord(2,:), '.r');

% xlim(x_lim);
% ylim(y_lim);

% spk_idx_L = nearest_idx3(S_L.t{1}, linpos.L.tvec);
% spk_idx_R = nearest_idx3(S_R.t{1}, linpos.R.tvec);
%
% subplot(2,2,4)
% c_ord = linspecer(6);
% cla
% hold on
% rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 0,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 25], 'facecolor', [c_ord(2,:), .2], 'edgecolor', 'none')
% rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 25,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 5], 'facecolor', [c_ord(6,:), .2], 'edgecolor', 'none')
%
% rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 30,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 20], 'facecolor', [c_ord(3,:), .2], 'edgecolor', 'none')
% rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 50,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 20], 'facecolor', [c_ord(5,:), .2], 'edgecolor', 'none')
% rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 70,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 20], 'facecolor', [c_ord(4,:), .2], 'edgecolor', 'none')
% rectangle('position', [min([linpos.L.tvec, linpos.R.tvec]), 90,  max([linpos.L.tvec, linpos.R.tvec]) - min([linpos.L.tvec, linpos.R.tvec]), 10], 'facecolor', [c_ord(1,:), .2], 'edgecolor', 'none')
%
% hlinl = plot(linpos.L.tvec, linpos.L.data, '.', 'color', [0.2 0.2 .2]);
% hlinls = plot(linpos.L.tvec(spk_idx_L), linpos.L.data(spk_idx_L), 'x', 'color', [0 0 1]);
%
% xlim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])])
% hlinr = plot(linpos.R.tvec, linpos.R.data, '.', 'color', [.2 0.2 0.2]);
% hlinrs = plot(linpos.R.tvec(spk_idx_R), linpos.R.data(spk_idx_R), 'x', 'color', [1 0 0]);
%
% xlim([min([linpos.L.tvec, linpos.R.tvec]), max([linpos.L.tvec, linpos.R.tvec])]);
% legend([hlinl, hlinr, hlinls, hlinrs], 'Linearized left trajectories', 'Linearized right trajectories', 'Left spikes', 'Right spikes');
% yyaxis right
% cla
% ylim([0 100])
% set(gca, 'ytick', [12.5 27.5 40 60 80 95], 'YTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})

%% generate TCs

Fs = 30;
occ_L = histc(linpos.L.data, bins)/Fs;
occ_L = trim_histc(occ_L);
occ_L(occ_L==0) = NaN;
occ_L = conv(occ_L, gkernal, 'same');

occ_R = histc(linpos.R.data, bins)/Fs;
occ_R = trim_histc(occ_R);
occ_R(occ_R==0) = NaN;
occ_R = conv(occ_R, gkernal, 'same');
%%
if cfg.plot
    figure(102)
    clf
end
for iS = 1:length(S.t)
    % overall spike stats
    fr = length(S.t{iS}) / (pos.tvec(end) - pos.tvec(1));
    binsize = 0.1; % select a small bin size for good time resolution
    tbin_edges = pos.tvec(1):binsize:pos.tvec(end);
    tbin_centers = tbin_edges(1:end-1)+binsize/2;
    
    
    spk_count = histc(S.t{iS},tbin_edges);
    spk_count = spk_count(1:end-1);
    gau_sdf = conv2(spk_count,gkernal,'same'); % convolve with gaussian window
    
    [~, fr_mean, fr_sd] = zscore(gau_sdf);
    
    
    %     plot([S.t{iS} S.t{iS}],[-1 -0.5],'Color',[0 0 0]); % note, plots all spikes in one command
    %     hold on;
    %     h = bar(tbin_centers,spk_count);
    %     set(h,'BarWidth',1,'EdgeColor','none','FaceColor',[1 0 0]); % reformat bar appearance
    %     plot(tbin_centers,gau_sdf,'g');
    %%
    
    if isempty(S_L.t{iS})
        TCs{iS}.spk_L = [];
        TCs{iS}.spk_tc_L = [];
        TCs{iS}.spk_tc_L_z = [];
        spk_idx_L = [];
    else
        % spk_hist
        TCs{iS}.spk_L = interp1(linpos.L.tvec, linpos.L.data, S_L.t{iS}, 'linear');
        TCs{iS}.spk_h_L = histc(TCs{iS}.spk_L', bins);
        TCs{iS}.spk_h_L = conv(TCs{iS}.spk_h_L(1:end-1), gkernal, 'same');
        TCs{iS}.spk_tc_L = TCs{iS}.spk_h_L./occ_L;
        
        TCs{iS}.spk_tc_L_z = (TCs{iS}.spk_tc_L - fr_mean) ./ fr_sd;
        spk_idx_L = nearest_idx3(S_L.t{iS}, linpos.L.tvec);
        
    end
    
    if isempty(S_R.t{iS})
        TCs{iS}.spk_R = [];
        TCs{iS}.spk_tc_R = [];
        TCs{iS}.spk_tc_R_z = [];
        spk_idx_R = [];
    else
        TCs{iS}.spk_R = interp1(linpos.R.tvec, linpos.R.data, S_R.t{iS}, 'linear');
        TCs{iS}.spk_h_R = histc(TCs{iS}.spk_R', bins);
        TCs{iS}.spk_h_R = conv(TCs{iS}.spk_h_R(1:end-1), gkernal, 'same');
        TCs{iS}.spk_tc_R = TCs{iS}.spk_h_R./occ_R;
        
        TCs{iS}.spk_tc_R_z = (TCs{iS}.spk_tc_R - fr_mean) ./ fr_sd;
        spk_idx_R = nearest_idx3(S_R.t{iS}, linpos.R.tvec);
    end
    %% plot TCs
    if cfg.plot
        
        m = length(S.t);
        n = 5;
        
        subplot(n,m,[iS m+iS m*2+iS])
        c_ord = linspecer(6);
        cla
        hold on
        if isempty(linpos.L.tvec)
            tvec_L = [];
        else
            tvec_L = linpos.L.tvec - min([linpos.L.tvec(1) linpos.R.tvec(1)]) ;
        end
        if isempty(linpos.R.tvec)
            tvec_R = [];
        else
            tvec_R = linpos.R.tvec - min([linpos.L.tvec(1) linpos.R.tvec(1)]) ;
        end

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
        % legend([hlinl, hlinr, hlinls, hlinrs], 'Linearized left trajectories', 'Linearized right trajectories', 'Left spikes', 'Right spikes', 'location', 'northoutside');
        if isempty(spk_idx_L)
            legend([ hlinrs], 'Right spikes', 'location', 'northoutside', 'Orientation', 'horizontal');
        elseif isempty(spk_idx_R)
            legend([ hlinls,], 'Left spikes', 'location', 'northoutside', 'Orientation', 'horizontal');
            if ~isempty(spk_idx_L) &&  ~isempty(spk_idx_R)
                legend([ hlinls, hlinrs], 'Left spikes', 'Right spikes', 'location', 'northoutside', 'Orientation', 'horizontal');
            end
        end
        set(gca, 'xtick', [12.5 27.5 40 60 80 95], 'xTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})
        xtickangle(25)
        if ~isempty(S_L.t{iS})
            subplot(n,m,[m*3+iS])
            cla
            hold on
            yyaxis left
            plot(bin_c, occ_L,  '-', 'color', [.7 .7 .7])
            ylim([min([occ_L occ_R]), max([occ_L occ_R])]);
            set(gca, 'YColor', 'k')
            ylabel('Occupancy L (s)')
            yyaxis right
            plot(bin_c,TCs{iS}.spk_tc_L, '-', 'color', c_ord(2,:))
            ylim([min([TCs{iS}.spk_tc_L TCs{iS}.spk_tc_R]), max([TCs{iS}.spk_tc_L TCs{iS}.spk_tc_R])]);
            ylabel('firing rate (Hz)')
            set(gca, 'YColor', c_ord(1,:))
            set(gca, 'xtick', [12.5 27.5 40 60 80 95], 'xTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})
            xtickangle(25)
        end
        
        if ~isempty(S_R.t{iS})
            subplot(n,m,[m*4+iS])
            cla
            hold on
            yyaxis left
            plot(bin_c, occ_R,  '-', 'color', [.7 .7 .7])
            ylim([min([occ_L occ_R]), max([occ_L occ_R])]);
            set(gca, 'YColor', 'k')
            ylabel('Occupancy R (s)')
            yyaxis right
            plot(bin_c,TCs{iS}.spk_tc_R, '-', 'color', c_ord(1,:))
            ylim([min([TCs{iS}.spk_tc_L TCs{iS}.spk_tc_R]), max([TCs{iS}.spk_tc_L TCs{iS}.spk_tc_R])]);
            ylabel('firing rate (Hz)')
            set(gca, 'YColor', c_ord(1,:))
            set(gca, 'xtick', [12.5 27.5 40 60 80 95], 'xTickLabel', {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'})
            xtickangle(25)
        end
    end
end
%% collect outputs

TC_out = [];
if exist('TCs')
    TC_out.TCs = TCs;
    TC_out.bin_c = bin_c;
    TC_out.bin_ticks_lines = [0 25 30 50 70 90 100]; 
    TC_out.bin_ticks =TC_out.bin_ticks_lines(1:end-1)+ (diff(TC_out.bin_ticks_lines)/2);  %[12.5 27.5 40 60 80 95];
    TC_out.bin_tick_labels = {'Start arm','Decision' 'Arm 1', 'Arm 2', 'Arm 3', 'Reward'};
end

