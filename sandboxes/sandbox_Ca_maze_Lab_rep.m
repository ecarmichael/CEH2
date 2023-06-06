%% sandbox_W maze checks

 inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\Raw\Inter';
% Load some data

load([inter_dir filesep 'dSub_g8_B_2023_02_03_maze.mat'])
figure(10)
plot(Maze.pos.data(1,:), Maze.pos.data(2,:), '.k')
Maze.pos.data(1,:) = Maze.pos.data(3,:)./7.8; 
Maze.pos.data(2,:) = Maze.pos.data(4,:)./8.9; 
Maze.pos.data = Maze.pos.data(1:2,:); 

Maze.pos.label = {'x', 'y'}; 

figure(10)
plot(Maze.pos.data(1,:), Maze.pos.data(2,:), '.k')
Maze.pos.data(1,:) = smoothdata(Maze.pos.data(1,:),'gaussian', round(1/mode(diff(Maze.pos.tvec)))*2); 
Maze.pos.data(2,:) = smoothdata(Maze.pos.data(2,:),'gaussian', round(1/mode(diff(Maze.pos.tvec)))*2); 
hold on
plot(Maze.pos.data(1,:), Maze.pos.data(2,:), '.r')


linspeed = getLinSpd([],Maze.pos); % linear speed 
% Threshold speed
cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = [5 ]; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg,linspeed); % only keep intervals with speed above thresh%% check plot

pos_move = restrict(Maze.pos, iv_fast); 

MS_Ca_check(Maze.ms)

c_ord = linspecer(8); 
%% activity by position


figure(101)
clf
ax(1) = subplot(7,1,1:5);
%     imagesc( Maze.ms.tvec,1:Maze.ms.numNeurons,  zscore(Maze.ms.denoise)'); set(gca, 'YDir', 'normal'); 
MS_Ca_Raster(Maze.ms.Binary', Maze.ms.tvec', 1); 
    title('Binarized traces')
    ylabel('cell ID')
     c_val = caxis; 
    caxis([0 16])
%     set(gca,'ytick', 0:100:ms.numNeurons*mult_fac, 'YTickLabel', (0:100:length(ms.units)*mult_fac)/mult_fac, 'TickDir', 'out')
    ylim([0 Maze.ms.numNeurons])
%     xlim([ms.time(1) ms.time(end)])

set(gca, 'xtick', []); 

ax(2) = subplot(7,1,6);
cla
hold on
plot(Maze.ms.tvec, mean(Maze.ms.denoise, 2), 'k')
gscatter(Maze.ms.tvec, mean(Maze.ms.denoise, 2), round(Maze.HD.data(1,:),2))
% colorbar('Location', 'northoutside', 'orientation' , 'horizontal')
ylabel({'Mean cell'; 'activity'})
legend('off')
hold on
y_lim = ylim; 
for ii = 1:length(maze.trials.types)
   if contains(maze.trials.types{ii},   'R correct')
       rectangle('position', [maze.trials.times(ii,1), y_lim(2)-0.005, diff(maze.trials.times(ii,:)) 0.005], 'FaceColor', c_ord(4,:));
   end
   if contains(maze.trials.types{ii},   'L correct')
       rectangle('position', [maze.trials.times(ii,1), y_lim(2)-0.005, diff(maze.trials.times(ii,:)) 0.005], 'FaceColor', c_ord(5,:));
   end
    
end
set(gca, 'xtick', []); 



% vline(maze.events.

ax(3) = subplot(7,1,7);
cla
yyaxis left
hold on
plot(Maze.pos.tvec, Maze.pos.data(1,:)+100,'-', 'color', c_ord(1,:))
plot(Maze.pos.tvec, Maze.pos.data(2,:)+100,'-', 'color', c_ord(2,:))
ylim([50 200])
ylabel('Position (cm)')
% set(gca, 'xtick', []); 


yyaxis right
plot(linspeed.tvec, linspeed.data(1,:))
ylim([0 50])
ylabel({'Velocity (cm/s)'})
%     xlabel('time (s)')
legend({'x', 'y', 'velocity'}, 'Orientation', 'horizontal')
% set(gca, 'xtick', []); 
% HD
% ax(4) = subplot(7,1,7);
% cla
% hold on
% plot(Maze.HD.tvec, Maze.HD.data(1,:), 'color', c_ord(5,:))
xlabel('time')
% 
% legend({'HD'})

linkaxes(ax, 'x')
xlim([Maze.ms.tvec(1) Maze.ms.tvec(end)])
% xlim([65 260])


%% linearize and trialify
cfg.bin_s = 2.5;
cfg.gk_win = 5;
cfg.gk_sd = 2.5;

bins = 0-(cfg.bin_s/2):cfg.bin_s:100+(cfg.bin_s/2);
bin_c = bins(1:end-1)+cfg.bin_s/2;
gkernal = gausskernel(cfg.gk_win/cfg.bin_s, cfg.gk_sd/cfg.bin_s);

cfg_c = [];
cfg_c.binsize = 1;
cfg_c.run_dist = 135;
coordM_cm = []; coordM_cm.coord = Common_CoorD.CoorD_N.coord; coordM_cm.units = 'cm';

coord.N = StandardizeCoord(cfg_c,coordM_cm, cfg_c.run_dist);


% %
% cfg_l = [];
% cfg_l.Coord = coord.L;
linpos.N = LinearizePos([],pos_move, coord.N);

Fs = ceil(1/mode(diff(Maze.pos.tvec)));
occ_M = histc(linpos.N.data, bins)/Fs;
occ_M = trim_histc(occ_M);
occ_M(occ_M==0) = NaN;
occ_M = conv(occ_M, gkernal, 'same');

%% prepare cells
    % grosmark method
    Csp = Maze.ms.deconv./Maze.ms.denoise; 
    Csp = Csp > 0.01; 
    Maze.ms.Csp = Csp; 
    
    binsize = 0.1; % in seconds, so everything else should be seconds too
    gauss_window = 1./binsize; % 1 second window
    gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
    gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
    gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window
    
    gau_z = zscore(gau_sdf, [], 2);


%% loop over cells
cfg_s.units = 'ms';
% S = MS_Binary2TS(cfg_s, Maze.ms);
S = MS_Deconv2TS(cfg_s, Maze.ms);
S_r = restrict(S, iv_fast)


binsize = 0.1; % select a small bin size for good time resolution
tbin_edges = Maze.pos.tvec(1):binsize:Maze.pos.tvec(end);
tbin_centers = tbin_edges(1:end-1)+binsize/2;
Q_mat = [];
for iS = size(gau_z,2):-1:1
    

    
    
%     spk_count = histc(S.t{iS},tbin_edges);
%     spk_count = spk_count(1:end-1);
%     gau_sdf = conv2(spk_count,gkernal,'same'); % convolve with gaussian window
%     
%     [~, fr_mean, fr_sd] = zscore(gau_sdf);

        TCs{iS}.spk_M = interp1(linpos.N.tvec, linpos.N.data, S.t{iS}, 'linear');
        TCs{iS}.spk_h_M = histc(TCs{iS}.spk_M', bins);
        Q_mat(iS,:) = histc(S.t{iS}, tbin_edges); 
        TCs{iS}.spk_h_M = conv(TCs{iS}.spk_h_M(1:end-1), gkernal, 'same');
        TCs{iS}.spk_tc_M = TCs{iS}.spk_h_M./occ_M;
  
        TC(iS,:) = TCs{iS}.spk_tc_M; 
        
        % 2d rate map
%         coord_vals = 1:length(coord.N.coord);
%         pos_grid = Maze.pos; 
%         pos_grid.data = interp2(coord.N.coord(1,:),coord.N.coord(2,:),coord_vals,Maze.pos.data(1,:),Maze.pos.data(2,:),'nearest')        
    
%         imagesc(coord.N.coord(1,:),coord.N.coord(2,:), occ_M)
end

% Q_

%% visualize Q
% [~, peak_idx] = max(Q_mat, [], 2);
s_bins = (1:100)*1.35; 

[~, l_idx] = max(TC,[],2); 
[~, ls_idx] = sort(l_idx); 
TC_s = TC(ls_idx,:); 
TC_m_s = TC_s./max(TC_s,[],2); 

figure(909)
clf
subplot(4,1,1:3)
imagesc(bin_c, 1:size(TC,2), TC_s)
% set(gca, 'YDir', 'normal')
% rectangle('position', [-2, 0, 51, 1], 'FaceColor', c_ord(4,:))
% rectangle('position', [26, 0, 3, 1], 'FaceColor', c_ord(1,:))
% rectangle('position', [49, 0, 6, 1], 'FaceColor', c_ord(2,:))
% rectangle('position', [55, 0, 46, 1], 'FaceColor', c_ord(3,:))
% rectangle('position', [77, 0, 3, 1], 'FaceColor', c_ord(5,:))
set(gca, 'xtick', [])
vline(23.5, 'w', 'R Reward')
vline(74, 'w', 'L Reward')
vline(0.6, 'w', 'R arm')
vline(53, 'w', 'L arm')
vline(46, 'w', 'C arm')
ylabel('Cell ID (peak TC sorted)')
xlim([0.5 (bin_c(end))+.5])

subplot(4,1,4)
histogram(l_idx)
xlim([0.5 length(bin_c)+.5])
% xline(29*1.35, 'color',c_ord(4,:))
% xline(31*1.35, 'color',c_ord(4,:))
% xline(49*1.35, 'color',c_ord(3,:))
% ylim([0 size(Q_mat,2)])
x_t = get(gca, 'xtick');
set(gca,'XTickLabel', round(bins(x_t)*1.35))
xlabel('distance on track (cm)')

%% 2d tuning
SET_xmin = 0; SET_ymin = 0; % set up bins
SET_xmax = 45; SET_ymax = 80;
SET_xBinSz = 2.5; SET_yBinSz = 2.5;

kernel = gausskernel([1 1],1); % 2d gaussian in bins

x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

occ_hist = histcn(pos_move.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');
 
no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
occ_hist(no_occ_idx) = NaN;
 
occ_hist = occ_hist .* (1/Fs); % convert samples to seconds using video frame rate (30 Hz)

n = 4;
m = 8;
s1_idx = 1:2:n*m; 
s2_idx = 2:2:n*m; 

ip = 0;
for iS = 1:size(gau_z,2)
    if length(S.t{ls_idx(iS)}) < 10
        continue
    end
spk_x = interp1(pos_move.tvec,pos_move.data(1,:),S.t{ls_idx(iS)},'linear');
spk_y = interp1(pos_move.tvec,pos_move.data(2,:),S.t{ls_idx(iS)},'linear');

    if ip >= (n*m)/2
        figure(102+iS)
        ip = 1;
    else
        ip = ip+1;
    end
    
    disp(iS)
    subplot(m, n, s1_idx(ip))
    hold on
    plot(pos_move.data(1,:), pos_move.data(2,:), '.','Color',[0.5 0.5 0.5],'MarkerSize',1)
    scatter(spk_x,spk_y,ones(1,length(spk_x))*20, S.t{ls_idx(iS)}, 'filled');
    axis off
    
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx) = NaN;
    tc = spk_hist./occ_hist;
    
subplot(m, n, s2_idx(ip))
pcolor(tc'); shading flat; axis off; %colorbar('Location', 'northoutside')

end

%%

Ass_Temp = assembly_patterns(gau_z(:,1:200));

time_proj = assembly_activity(Ass_Temp,gau_sdf(:,1:200)'); 


%% color code the assemblies 
Ass_sort = []; 
for ii = size(Ass_Temp,2):-1:1
    
    [Ass_sort(:,ii), this_idx] = sort(Ass_Temp(:,ii), 'descend'); 
    thresh = prctile(Ass_Temp(:,ii), 95); 
    c_idx = Ass_Temp(this_idx,ii) > thresh; 
    
%     stem(Ass_Temp(c_idx, ii))
    Ass_col(ii,:) = c_idx;
    
end

%% stem plot for first few ensembles 
 figure(303);clf;  hold on
for ii = 1:5
        keep_idx  = Ass_Temp(:,ii) > .2; 
        
        stem(Ass_Temp(:,ii)); 
%                 stem(Ass_Temp(keep_idx,ii))

end


%%
figure(202)
ax(1) = subplot(6, 1, 1)
hold on
plot(Maze.ms.tvec, Maze.pos.data(1,:)); 
plot(Maze.ms.tvec, Maze.pos.data(2,:)); 
xlim([Maze.ms.tvec(1) Maze.ms.tvec(end)])
ylabel('position')
    set(gca, 'xtick', []);

ax(2) = subplot(6,1,2:5);
hold on
for ii = size(gau_sdf(:,1:200), 2):-1:1
    s_t = Maze.ms.tvec(Csp(:,ii) >0);
    if ~isempty(s_t)
        if ii == 1
           plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', 'k', 'linewidth', 4)
           plot([s_t, s_t]', [(ones(size(s_t())*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', 'k', 'linewidth', 4)

        else
            
            plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', 'k', 'linewidth', 4)
        end
    end
%     plot(ms_trk.time, gau_sdf(:,ii)+ii)
end
ylabel('Cell ID')
xlim([Maze.ms.tvec(1) Maze.ms.tvec(end)])
ylim([0 size(Csp(:,1:200), 2)])
        set(gca, 'xtick', []);

ax(3) = subplot(6,1,6);
cla
hold on
c_ord = linspecer(5); 
for ii = 1:3
   plot3(Maze.ms.tvec, zscore(time_proj(ii,:)),ones(1, length(time_proj))*100* ii, 'color', c_ord(ii,:)) 
    
end
    xlim([Maze.ms.tvec(1) Maze.ms.tvec(end)])
ylabel({'Reactivation'; 'strength'})
xlabel('time (s)')
linkaxes(ax, 'x')
legend({'1', '2', '3'}, 'Orientation', 'horizontal', 'Location', 'northeast', 'Box', 'off')
