function sandbox_v4_dSub_space_screener(cfg_in, fname)
%% MS_dSub_space_screener:
%
%
%
%    Inputs:
%    - cfg [struct]   configuration see the defaults below.
%
%    - data_dir:  [string]   path to nlx data folder. default is current
%    directory if empty.
%
%    - save_dir: [string] path to save output plots and files. defaults to
%    same as data_dir if empty
%
%    Outputs:
%    -
%
%
%
%
% EC 2022-10-03   initial version
%
%
%
%% initialize
% if nargin == 0
%     cfg_in = [];
%     data_dir = cd;
%     save_dir = data_dir;
% elseif nargin == 1
%     data_dir = cd;
%     save_dir = data_dir;
% elseif nargin == 2
%     save_dir = data_dir;
% end
% nice colours
c_ord = linspecer(4);
cfg_def = [];



cfg = ProcessConfig(cfg_def, cfg_in);
%%  Load and split recordings
parts = strsplit(fname, '_');
type = parts{end}; 
sess = [parts{end-3} '_' parts{end-2} '_' parts{end-1}]; 
subject =  strsplit(fname, '_2023');
subject = subject{1}; 


% load('maze.mat');

% load([fileparts(data_dir) filesep 'Common_CoorD.mat'])
% 
% if str2double(sess(2:end)) < 14
%     Common_CoorD.CoorD_L.coord(1,:) = Common_CoorD.CoorD_L.coord(1,:) - 4;
%     Common_CoorD.CoorD_R.coord(1,:) = Common_CoorD.CoorD_R.coord(1,:) - 4;
% end
%% load the maze data
load([fname '.mat'])
if strcmpi(type, 'maze')
    S = Maze.minianms.S;
    ms = Maze.minianms;
    pos_w = Maze.pos;
end
pos_w.data = pos_w.data(1:2,:); 
pos_w.label = [];
pos_w.label = {'x', 'y'}; 

linspeed = getLinSpd([],pos_w); % linear speed

% Threshold speed
cfg_l = []; cfg_l.method = 'raw'; cfg_l.operation = '>'; cfg_l.threshold = 2.5; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg_l,linspeed); % only keep intervals with speed above thresh

% Restrict data so it includes fast intervals only
pos_w = restrict(pos_w,iv_fast);
S_w = restrict(S,iv_fast);

%% if it exists load the open field

if exist(strrep([fname '.mat'], 'maze', 'OF'),'file')
load(strrep([fname '.mat'], 'maze', 'OF'))

S = OF.minianms.S; 
pos_o = OF.pos;
pos_o.data = pos_o.data(1:2,:); 
pos_o.label = [];
pos_o.label = {'x', 'y'}; 

linspeed = getLinSpd([],pos_o); % linear speed

% Threshold speed
cfg_l = []; cfg_l.method = 'raw'; cfg_l.operation = '>'; cfg_l.threshold = 2.5; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg_l,linspeed); % only keep intervals with speed above thresh

% Restrict data so it includes fast intervals only
pos_o = restrict(pos_o,iv_fast);
S_o = restrict(S,iv_fast);
end


%% Spatial tuning in the maze
% % trial times
% R_idx = contains( maze.trials.types, {'FR correct'; 'CR correct'; 'CL error'});
% L_idx = contains( maze.trials.types, {'FL correct'; 'CL correct'; 'CR error'});
% 
% 
% if length(maze.events.Box_in(:,1)) < length(maze.trials.times(:,2))
%     % add box to trials
%     R_idx = R_idx(1:length(maze.events.Box_in(:,1)));
%     L_idx = L_idx(1:length(maze.events.Box_in(:,1)));
% 
%     W_iv = iv(maze.events.Box_in(:,1), maze.trials.times(1:length(maze.events.Box_in(:,1)),2)+5);
%     R_iv = iv(maze.events.Box_in(R_idx,1), maze.trials.times(R_idx,2)+5);
%     L_iv = iv(maze.events.Box_in(L_idx,1), maze.trials.times(L_idx,2)+5);
%     
% elseif length(maze.events.Box_in(:,1)) > length(maze.trials.times(:,2))
%      R_idx = R_idx(1:length(maze.trials.times(:,2)));
%     L_idx = L_idx(1:length(maze.trials.times(:,2)));
%     
%     W_iv = iv(maze.events.Box_in(1:length(maze.trials.times(:,1)),1), maze.trials.times(:,2)+5);
%     R_iv = iv(maze.events.Box_in(R_idx,1), maze.trials.times(R_idx,2)+5);
%     L_iv = iv(maze.events.Box_in(L_idx,1), maze.trials.times(L_idx,2)+5);
% else
%     
%     % add box to trials
%     W_iv = iv(maze.events.Box_in(:,1), maze.trials.times(:,2)+5);
%     R_iv = iv(maze.events.Box_in(R_idx,1), maze.trials.times(R_idx,2)+5);
%     L_iv = iv(maze.events.Box_in(L_idx,1), maze.trials.times(L_idx,2)+5);
% end
% 
% 
% 
% pos_W = restrict(pos_w, W_iv);
% pos_L = restrict(pos_w, L_iv);
% pos_R = restrict(pos_w, R_iv);
% 
% S_W = restrict(S_w, W_iv);
% S_L = restrict(S_w, L_iv);
% S_R = restrict(S_w, R_iv);
% 
% %% grab the linearized maze data TCs
% cfg_tc = [];
% cfg_tc.plot = 0;
% TC = dSub_gen_1d_TC(cfg_tc, cd, Common_CoorD);

%% generate rate map

for iS = 1:length(S.t)
    %%
    figure(iS)
    clf
    
    % get the spiking
    spk_x = interp1(pos_w.tvec,pos_w.data(1,:),S_w.t{iS},'linear');
    spk_y = interp1(pos_w.tvec,pos_w.data(2,:),S_w.t{iS},'linear');
    spk_idx = nearest_idx3(S_w.t{iS}, pos_w.tvec); 
    tvec_cord = winter(length(pos_w.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
    
    ax1= subplot(4,6,[1:3 7:9]);
    scatter(pos_w.data(1,:), pos_w.data(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    colormap(ax1, 'winter');
    cax = colorbar;
    %     cax.Position(1) = cax.Position(1) + 0.03;
    cax.Ticks = [0 1]; cax.TickLabels = {'start', 'end'};
    
    hold on
    hr =  plot(spk_x,spk_y, '.r', 'markersize', 6);
    
%     plot(Common_CoorD.CoorD_L.coord(1,:), Common_CoorD.CoorD_L.coord(2,:), '--', 'color', c_ord(1,:))
%     plot(Common_CoorD.CoorD_R.coord(1,:), Common_CoorD.CoorD_R.coord(2,:), '--', 'color', c_ord(2,:))
    ax2= subplot(4,6,[13:15]);
    hold on
    scatter(pos_w.tvec, pos_w.data(1,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    scatter(pos_w.tvec, pos_w.data(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    plot(pos_w.tvec(spk_idx),pos_w.data(1,spk_idx) , '.r', 'markersize', 6);
    plot(pos_w.tvec(spk_idx),pos_w.data(2,spk_idx) , '.r', 'markersize', 6);
    
    % rate map
    if max(pos_w.data(1,:)) <300 && max(pos_w.data(2,:)) < 300
        SET_xmin = 25; SET_ymin = 25; % set up bins
        SET_xmax = 275; SET_ymax = 275;
        SET_xBinSz = 10; SET_yBinSz =10;
    else
        SET_xmin = 50; SET_ymin = 50; % set up bins
        SET_xmax = 600; SET_ymax = 400;
        SET_xBinSz = 10; SET_yBinSz =10;
    end
    
%     NaN_map = zeros(length(SET_xmin:SET_xBinSz:SET_xmax), length(SET_ymin:SET_yBinSz:SET_ymax));
%     NaN_map(1:5,25:end) = NaN;
%     NaN_map(30:end,25:end) = NaN;
%     NaN_map(:,29:end) = NaN;
%     nan_idx = isnan(NaN_map);
    
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    % set up gaussian
    kernel = gausskernel([SET_xBinSz SET_yBinSz],SET_xBinSz/2); % 2d gaussian in bins
    
    
    % compute occupancy
    occ_hist = hist3(pos_w.data(1:2,:)', 'edges', {x_edges y_edges});
    %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    
    no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
    occ_hist = conv2(occ_hist,kernel,'same');
    occ_hist(no_occ_idx) = NaN;
    %      occ_hist(nan_idx) = NaN;
    
    occ_hist = occ_hist .* mode(diff(pos_w.tvec)); % convert samples to seconds using video frame rate (30 Hz)
    
    %     subplot(2,2,3)
    %     pcolor(occ_hist'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'secs'; cb.Ticks = [0 cb.Ticks(end)];
    %     title('occupancy');
    
    % get the spike map
    spk_hist = hist3([spk_x, spk_y], 'edges', {x_edges y_edges});
    %     spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    
    spk_hist = conv2(spk_hist,kernel,'same');
    spk_hist(no_occ_idx) = NaN;
    
    %     subplot(2,2,3)
    %     pcolor(spk_hist'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'nSpikes'; cb.Ticks = [0 cb.Ticks(end)];
    %     title('spikes');
    
    % rate map
    tc = spk_hist./occ_hist;
    
    subplot(4,6,[4:6 10:12]);
    pcolor(tc'); shading flat; axis off;
    cb=colorbar; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)];
    
    % plot the cell info
    subplot(4,6,[17])   
    imagesc(mean(ms.SFP(:,:,:),3))
    hold on

        [row, col] = find(ms.SFP(:,:,iS) == max(max(ms.SFP(:,:,iS))));
        scatter(col, row,50, 'MarkerEdgeColor', c_ord(2,:), 'LineWidth', 2)
    title(['nCells: ' num2str(length(ms.units(iS)))]);
    xlim([min(ms.centroids(1,:))-min(ms.centroids(1,:))*.2  max(ms.centroids(1,:))+max(ms.centroids(1,:))*.2])
    ylim([min(ms.centroids(2,:))-min(ms.centroids(2,:))*.2  max(ms.centroids(2,:))+max(ms.centroids(2,:))*.2])

        subplot(4,6,[22 23 24])   
        cla
        hold on
                plot(ms.tvec, ms.Deconv(:,iS), 'color', c_ord(1,:));

        plot(ms.tvec, ms.RawTraces(:,iS), 'color', [.8 .8 .8 .2]); 
xlim([ms.tvec(1), ms.tvec(end)])
    
%     
%     % plot the linear TCs
%     if ~isempty(TC)
%         subplot(2,6,7:8);
%         this_tc = TC.TCs{iS};
%         if isempty(this_tc.spk_tc_L)
%             this_tc.spk_tc_L = nan(size(this_tc.spk_tc_R)); 
%         end
%         if isempty(this_tc.spk_tc_R)
%             this_tc.spk_tc_R = nan(size(this_tc.spk_tc_L)); 
%         end
%         imagescnan(TC.bin_c, 1:2, [this_tc.spk_tc_L; this_tc.spk_tc_R])
%         cb = colorbar;
%         set(gca, 'XTick', TC.bin_ticks, 'XTickLabel', TC.bin_tick_labels, 'ytick', 1:2, 'YTickLabel', {'Right', 'Left'}, 'color', 'k');
%         vline(TC.bin_ticks_lines(2:end), 'r'); 
%         cb=colorbar; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)];
%         xtickangle(45)
%     end
%     %% OF rate map
%     
%     % interpolate tracking off the box floor.
%     OF_pos = pos_OF.data;
%     x_rm_idx = OF_pos(1,:) > 85 | OF_pos(1,:) < 42;
%     y_rm_idx = OF_pos(2,:) > 75 | OF_pos(2,:) < 38;
%     OF_pos(1,x_rm_idx) = NaN;
%     OF_pos(2,y_rm_idx) = NaN;
%     
%     OF_pos_int(1,:) = OF_pos(1,:);%fillmissing(OF_pos(1,:), 'linear');
%     OF_pos_int(2,:) = OF_pos(2,:); %fillmissing(OF_pos(2,:), 'linear');
%     
%     
%     %     figure(99)
%     %     clf
%     %     hold on
%     %     plot(pos_OF.data(1,:), pos_OF.data(2,:), '.k')
%     %     plot(OF_pos(1,:), OF_pos(2,:), '.r')
%     %     plot(OF_pos_int(1,:), OF_pos_int(2,:), '.b')
%     
%     
%     
%     spk_x = interp1(pos_OF.tvec,OF_pos_int(1,:),S_OF.t{iS},'linear');
%     spk_y = interp1(pos_OF.tvec,OF_pos_int(2,:),S_OF.t{iS},'linear');
%     
%     SET_xmin = 40; SET_ymin = 35; % set up bins
%     SET_xmax = 85; SET_ymax = 75;
%     SET_xBinSz = 5; SET_yBinSz =5;
%     
%     x_edges = SET_xmin:SET_xBinSz:SET_xmax;
%     y_edges = SET_ymin:SET_yBinSz:SET_ymax;
%     
%     % set up gaussian
%     kernel = gausskernel([SET_xBinSz SET_yBinSz],SET_xBinSz/2); % 2d gaussian in bins
%     
%     
%     % compute occupancy
%     occ_hist = hist3(OF_pos_int(1:2,:)', 'edges', {x_edges y_edges});
%     %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
%     
%     no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
%     occ_hist = conv2(occ_hist,kernel,'same');
%     occ_hist(no_occ_idx) = NaN;
%     %      occ_hist(nan_idx) = NaN;
%     
%     occ_hist = occ_hist .* mode(diff(pos_OF.tvec)); % convert samples to seconds using video frame rate (30 Hz)
%     
%     %     subplot(2,2,3)
%     %     pcolor(occ_hist'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'secs'; cb.Ticks = [0 cb.Ticks(end)];
%     %     title('occupancy');
%     
%     % get the spike map
%     spk_hist = hist3([spk_x, spk_y], 'edges', {x_edges y_edges});
%     %     spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
%     
%     spk_hist = conv2(spk_hist,kernel,'same');
%     spk_hist(no_occ_idx) = NaN;
%     
%     %     subplot(2,2,3)
%     %     pcolor(spk_hist'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'nSpikes'; cb.Ticks = [0 cb.Ticks(end)];
%     %     title('spikes');
%     
%     % rate map
%     tc = spk_hist./occ_hist;
%     
%     % running time as a colour
%     tvec_cord = winter(length(pos_OF.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
%     
%     ax2= subplot(2,6,[9 10]);
%     cla
%     scatter(OF_pos_int(1,:), OF_pos_int(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%     colormap(ax2, 'winter');
%     cax = colorbar;
%     %     cax.Position(1) = cax.Position(1) + 0.03;
%     cax.Ticks = [0 1]; cax.TickLabels = {'start', 'end'};
%     hold on
%     hr =  plot(spk_x,spk_y, '.r', 'markersize', 6);
%     axis off
%     
%     ax3= subplot(2,6,[11 12]);
%     pcolor(tc'); shading flat; axis off ; cb=colorbar; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)]; % cb.Position(1) = cb.Position(1)
%     
    %% print and save
%     if length(S_W.t{iS}) > 150
%         cfg_fig = [];
%         cfg_fig.ft_size = 12;
%         SetFigure(cfg_fig, gcf)
%         maximize
%         saveas(gcf, [save_dir filesep fname '_' S.label{iS}(1:end-2) '.png'])
%         % saveas(gcf, [save_dir filesep fname '_' S.label{iS}(1:end-2) '.fig'])
%         
%         this_TC = TC;
%         this_TC.TCs = this_TC.TCs{iS}; 
%         
%         save([save_dir filesep fname '_' S.label{iS}(1:end-2) '_tc.mat'], 'this_TC');
%     end
end
