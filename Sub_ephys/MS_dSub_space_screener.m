function MS_dSub_space_screener(cfg_in, data_dir, save_dir)
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
if nargin == 0
    cfg_in = [];
    data_dir = cd;
    save_dir = data_dir;
elseif nargin == 1
    data_dir = cd;
    save_dir = data_dir;
elseif nargin == 2
    save_dir = data_dir;
end

cfg_def = [];



cfg = ProcessConfig(cfg_def, cfg_in);
%%  Load and split recordings
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


load([fileparts(data_dir) filesep 'Common_CoorD.mat'])

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





pos_w = restrict(pos, blocks.W_maze(1), blocks.W_maze(2));
S_w = restrict(S, blocks.W_maze(1), blocks.W_maze(2));
evts_w = restrict(evts, blocks.W_maze(1), blocks.W_maze(2));

linspeed = getLinSpd([],pos_w); % linear speed

% Threshold speed
cfg_l = []; cfg_l.method = 'raw'; cfg_l.operation = '>'; cfg_l.threshold = 1; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg_l,linspeed); % only keep intervals with speed above thresh

% Restrict data so it includes fast intervals only
pos_w = restrict(pos_w,iv_fast);
S_w = restrict(S_w,iv_fast);


pos_OF = restrict(pos, blocks.OF(1), blocks.OF(2));
pos_OF.data = pos_OF.data*1.5;
S_OF = restrict(S, blocks.OF(1), blocks.OF(2));
evts_OF = restrict(evts, blocks.OF(1), blocks.OF(2));

% Threshold speed
linspeed = getLinSpd([],pos_OF); % linear speed
cfg_l = []; cfg_l.method = 'raw'; cfg_l.operation = '>'; cfg_l.threshold = 1; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg_l,linspeed); % only keep intervals with speed above thresh

% Restrict data so it includes fast intervals only
pos_OF = restrict(pos_OF,iv_fast);
S_OF = restrict(S_OF,iv_fast);



% nice colours
c_ord = linspecer(4);
%% Spatial tuning in the maze
% trial times
R_idx = contains( maze.trials.types, {'FR correct'; 'CR correct'; 'CL error'});
L_idx = contains( maze.trials.types, {'FL correct'; 'CL correct'; 'CR error'});

if length(maze.events.Box_in(:,1)) < length(maze.trials.times(:,2))
    % add box to trials
    R_idx = R_idx(1:length(maze.events.Box_in(:,1)));
    L_idx = L_idx(1:length(maze.events.Box_in(:,1)));

    W_iv = iv(maze.events.Box_in(:,1), maze.trials.times(1:length(maze.events.Box_in(:,1)),2)+5);
    R_iv = iv(maze.events.Box_in(R_idx,1), maze.trials.times(R_idx,2)+5);
    L_iv = iv(maze.events.Box_in(L_idx,1), maze.trials.times(L_idx,2)+5);
    
else
    
    % add box to trials
    W_iv = iv(maze.events.Box_in(:,1), maze.trials.times(:,2)+5);
    R_iv = iv(maze.events.Box_in(R_idx,1), maze.trials.times(R_idx,2)+5);
    L_iv = iv(maze.events.Box_in(L_idx,1), maze.trials.times(L_idx,2)+5);
end



pos_W = restrict(pos_w, W_iv);
pos_L = restrict(pos_w, L_iv);
pos_R = restrict(pos_w, R_iv);

S_W = restrict(S_w, W_iv);
S_L = restrict(S_w, L_iv);
S_R = restrict(S_w, R_iv);

%% grab the linearized maze data TCs
cfg_tc = [];
cfg_tc.plot = 0;
TC = dSub_gen_1d_TC(cfg_tc, cd, Common_CoorD);

%% generate rate map

for iS = 1:length(S.t)
    %%
    figure(iS)
    clf
    
    % get the spiking
    spk_x = interp1(pos_W.tvec,pos_W.data(1,:),S_W.t{iS},'linear');
    spk_y = interp1(pos_W.tvec,pos_W.data(2,:),S_W.t{iS},'linear');
    
    tvec_cord = winter(length(pos_W.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
    
    ax1= subplot(2,6,1:3);
    scatter(pos_W.data(1,:), pos_W.data(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    colormap(ax1, 'winter');
    cax = colorbar;
    %     cax.Position(1) = cax.Position(1) + 0.03;
    cax.Ticks = [0 1]; cax.TickLabels = {'start', 'end'};
    
    hold on
    hr =  plot(spk_x,spk_y, '.r', 'markersize', 6);
    
    plot(Common_CoorD.CoorD_L.coord(1,:), Common_CoorD.CoorD_L.coord(2,:), '--', 'color', c_ord(1,:))
    plot(Common_CoorD.CoorD_R.coord(1,:), Common_CoorD.CoorD_R.coord(2,:), '--', 'color', c_ord(2,:))
    % spike characteristics
    try % see if the waveform file exists and works.
        load([S.label{iS}(1:end-2) '-wv.mat'], 'mWV', 'xrange')
        
        subplot(2,6,1:3);
        hold on
        for ii  = 1:4
            plot ((xrange(:,ii)./5)- (xrange(1,1)./5),(mWV(:,ii)./2000)+50,"LineWidth",3, 'color', c_ord(ii,:));
        end
        %             legend("Ch 1", "Ch 2","Ch 3","Ch 4", 'fontsize', 8, 'location', 'north', 'orientation', 'horizontal' ); legend boxoff
        %         axis off
        title ({strrep(fname, '_', ' ') ;  strrep(S.label{iS}, '_', ' '); ['grade: ' S.label{iS}(end-2) ' (' num2str(length(S_W.t{iS})) ' spikes)'];})
        %         xlim([min(xrange,[], 'all') max(xrange,[], 'all')])
    catch
        %         subplot(2,2,1)
        %         text(0, .5, {'waveform file does not exist' ; 'or is corrupted'});
        %         axis off
    end
    legend([hr], {'spikes'})
    
    % rate map
    SET_xmin = 10; SET_ymin = 0; % set up bins
    SET_xmax = 100; SET_ymax = 70;
    SET_xBinSz = 5; SET_yBinSz =5;
    
    NaN_map = zeros(length(SET_xmin:SET_xBinSz:SET_xmax), length(SET_ymin:SET_yBinSz:SET_ymax));
    NaN_map(1:5,25:end) = NaN;
    NaN_map(30:end,25:end) = NaN;
    NaN_map(:,29:end) = NaN;
    nan_idx = isnan(NaN_map);
    
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    % set up gaussian
    kernel = gausskernel([SET_xBinSz SET_yBinSz],SET_xBinSz/2); % 2d gaussian in bins
    
    
    % compute occupancy
    occ_hist = hist3(pos_W.data(1:2,:)', 'edges', {x_edges y_edges});
    %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    
    no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
    occ_hist = conv2(occ_hist,kernel,'same');
    occ_hist(no_occ_idx) = NaN;
    %      occ_hist(nan_idx) = NaN;
    
    occ_hist = occ_hist .* mode(diff(pos_W.tvec)); % convert samples to seconds using video frame rate (30 Hz)
    
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
    
    subplot(2,6,4:6);
    pcolor(tc'); shading flat; axis off;
    cb=colorbar; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)];
    
    
    % plot the linear TCs
    if ~isempty(TC)
        subplot(2,6,7:8);
        this_tc = TC.TCs{iS};
        
        imagescnan(TC.bin_c, 1:2, [this_tc.spk_tc_L; this_tc.spk_tc_R])
        cb = colorbar;
        set(gca, 'XTick', TC.bin_ticks, 'XTickLabel', TC.bin_tick_labels)
        vline(TC.bin_ticks_lines(2:end), 'r'); 
        cb=colorbar; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)];
        xtickangle(45)
    end
    %% OF rate map
    
    % interpolate tracking off the box floor.
    OF_pos = pos_OF.data;
    x_rm_idx = OF_pos(1,:) > 85 | OF_pos(1,:) < 42;
    y_rm_idx = OF_pos(2,:) > 75 | OF_pos(2,:) < 38;
    OF_pos(1,x_rm_idx) = NaN;
    OF_pos(2,y_rm_idx) = NaN;
    
    OF_pos_int(1,:) = OF_pos(1,:);%fillmissing(OF_pos(1,:), 'linear');
    OF_pos_int(2,:) = OF_pos(2,:); %fillmissing(OF_pos(2,:), 'linear');
    
    
    %     figure(99)
    %     clf
    %     hold on
    %     plot(pos_OF.data(1,:), pos_OF.data(2,:), '.k')
    %     plot(OF_pos(1,:), OF_pos(2,:), '.r')
    %     plot(OF_pos_int(1,:), OF_pos_int(2,:), '.b')
    
    
    
    spk_x = interp1(pos_OF.tvec,OF_pos_int(1,:),S_OF.t{iS},'linear');
    spk_y = interp1(pos_OF.tvec,OF_pos_int(2,:),S_OF.t{iS},'linear');
    
    SET_xmin = 40; SET_ymin = 35; % set up bins
    SET_xmax = 85; SET_ymax = 75;
    SET_xBinSz = 5; SET_yBinSz =5;
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    % set up gaussian
    kernel = gausskernel([SET_xBinSz SET_yBinSz],SET_xBinSz/2); % 2d gaussian in bins
    
    
    % compute occupancy
    occ_hist = hist3(OF_pos_int(1:2,:)', 'edges', {x_edges y_edges});
    %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    
    no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
    occ_hist = conv2(occ_hist,kernel,'same');
    occ_hist(no_occ_idx) = NaN;
    %      occ_hist(nan_idx) = NaN;
    
    occ_hist = occ_hist .* mode(diff(pos_OF.tvec)); % convert samples to seconds using video frame rate (30 Hz)
    
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
    
    % running time as a colour
    tvec_cord = winter(length(pos_OF.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
    
    ax2= subplot(2,6,[9 10]);
    cla
    scatter(OF_pos_int(1,:), OF_pos_int(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    colormap(ax2, 'winter');
    cax = colorbar;
    %     cax.Position(1) = cax.Position(1) + 0.03;
    cax.Ticks = [0 1]; cax.TickLabels = {'start', 'end'};
    hold on
    hr =  plot(spk_x,spk_y, '.r', 'markersize', 6);
    axis off
    
    ax3= subplot(2,6,[11 12]);
    pcolor(tc'); shading flat; axis off ; cb=colorbar; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)]; % cb.Position(1) = cb.Position(1)
    
    %% print and save
    if length(S_W.t{iS}) > 150
        cfg_fig = [];
        cfg_fig.ft_size = 12;
        SetFigure(cfg_fig, gcf)
        maximize
        saveas(gcf, [save_dir filesep fname '_' S.label{iS}(1:end-2) '.png'])
        % saveas(gcf, [save_dir filesep fname '_' S.label{iS}(1:end-2) '.fig'])
        
        this_TC = TC;
        this_TC.TCs = this_TC.TCs{iS}; 
        
        save([save_dir filesep fname '_' S.label{iS}(1:end-2) '_tc.mat'], 'this_TC');
    end
end
