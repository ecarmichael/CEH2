function dsubscreener_ec_W_maze(data_dir, inter_dir)
%% dsubscreener:
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
% EP 2021-07-08   initial version
%
%
%
%% initialize
if nargin < 1
    data_dir = cd;
end
if nargin < 2
    inter_dir = [data_dir filesep 'inter'];
end

if ~exist(inter_dir, 'dir')
    mkdir(inter_dir)
end

c_ord = linspecer(5); % nice colours.

%% get session wide variables

% load the meta data
if isempty(dir('*meta.m'))
    MS_Write_meta_dSub;
end
Meta = MS_Load_meta;

% get position
cfg_pos = [];
cfg_pos.convFact = [8,8];
pos = MS_LoadPos(cfg_pos);
spd = getLinSpd([],pos);
% spd_gk = gausskernel(ceil(1/mode(diff(spd.tvec))),ceil(1/mode(diff(spd.tvec)))*.2);
% spd.data =  conv2(spd.data,gausswin(100),'same'); % smooth over 2 seconds
spd.data = smooth(spd.data, 2*ceil(1/mode(diff(spd.tvec))))';
spd.data(1) = spd.data(2); % correct for jump in first sample 

% get LFP
cfg_lfp = [];
cfg_lfp.fc = {Meta.goodCSC};
cfg_lfp.desired_sampling_frequency = 2000;
csc = MS_LoadCSC(cfg_lfp);


% get some file info. (can be done with Meta_keys later.
% assumes code is saved as 'Subject_yyyy_mm_dd_task'
fname = strsplit(data_dir, filesep);
fname = strrep(fname{end}, '-', '_');

% laod and split events
evt = LoadEvents([]);
start_idx = find(contains(evt.label, 'Starting Recording'));
end_idx = find(contains(evt.label, 'Stopping Recording'));




if length(evt.t{start_idx}) ~= 4
    rec_len = evt.t{end_idx} - evt.t{start_idx};
    [~, idx] = sort(rec_len, 'descend'); 
    keep_idx = sort(idx(1:4)); 
    
%     error('too many start and stop times. should be 4.  Figure this out')
else
    keep_idx = 1:4; 
end

block_idx = [];
if length(evt.t{start_idx}) > 1
    disp('splitting into blocks')
    Blocks = {'Pre_sleep', 'W_maze', 'OF', 'Post_sleep'};
    for iB = 1:length(Blocks)
        if strcmp(Blocks{iB}, 'Pre_sleep')  || strcmp(Blocks{iB}, 'Post_sleep') % trim sleep blocks. 
            block_idx.(Blocks{iB}) = [evt.t{start_idx}(keep_idx(iB)) evt.t{start_idx}(keep_idx(iB))+120*60];
        else
        block_idx.(Blocks{iB}) = [evt.t{start_idx}(keep_idx(iB)) evt.t{end_idx}(keep_idx(iB))];
        end
        fprintf('Block: %s duration = %.0fmins\n', Blocks{iB}, (block_idx.(Blocks{iB})(2) -block_idx.(Blocks{iB})(1))/60); 
    end
end

%% get all .t files

file_list = FindFiles('*.t');

for iC = 1:length(file_list)
    
    cfg = [];
    cfg.fc = file_list(iC);
    cfg.getTTnumbers = 0;
    
    S = LoadSpikes(cfg);
    
    % split into recording blocks.
    
    for iB = 1:length(Blocks)
        this_csc = restrict(csc, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
        this_S = restrict(S, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
        this_pos = restrict(pos, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
        this_spd = restrict(spd, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
%         this_spd.data(1) = this_spd.data(2); 
        this_evt = restrict(evt, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2)); 
        
        if isempty(this_S.t{1})
            continue
        end
        if iB == 2
            
%             trials = dSub_wmaze_trialfun([], this_pos)
            arms = {'l', 'c', 'r'}; maze.arms = arms; 
            for ii = 1:3
            if sum(strcmp(this_evt.label, arms{ii})) > 0
                maze.(arms{ii}) = this_evt.t{strcmp(this_evt.label, arms{ii})};
            else
                maze.(arms{ii}) = []; 
            end
            end
            fprintf('Left: %d   |  Center: %d   | Right: %d\n', length(maze.l), length(maze.c), length(maze.r))
            
        elseif iB == 3
                this_pos.data(1,:) = this_pos.data(1,:)/2;
                this_pos.data(2,:) = this_pos.data(2,:)/2;
                this_spd.data(1,:) = this_spd.data(1,:)*2;

        end
        
        % generate psd
        cfg_psd.hann_win = 2^11;
        if iB == 1 || iB ==4
            [Px, Fx] = pwelch(this_csc.data(1,floor(length(this_csc.data)/10):end), hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);
        else
            [Px, Fx] = pwelch(this_csc.data(1,:), hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);
        end
        figure((100*iC) + iB)
           SetFigure([], gcf);
        set(gcf, 'position', [81 -377  1760  880])
        subplot(3,4,9)
        hold on
        plot(Fx, 10*log10(Px), 'color', [244,93,1]/255, 'linewidth',2)
        xlim([0 150])
        set(gca,'XTick',0:20:140)
        xlabel("Frequency (Hz)")
        
        % add rectangles marking gamma and theta ranges
        color.blue = double([158,202,225])/255;
        color.green = double([168,221,181])/255;
        color.red = [0.9    0.3    0.3];
        
        y_val = ylim;
        
        % colour bars for specific frequencies of interest.
        rectangle('position', [6, y_val(1), 6, y_val(2) - y_val(1)],  'facecolor', [color.blue 0.3], 'edgecolor', [color.blue 0.3])
        rectangle('position', [30, y_val(1), 20, y_val(2) - y_val(1)],  'facecolor', [color.green 0.3], 'edgecolor', [color.green 0.3])
        rectangle('position', [70, y_val(1), 20, y_val(2) - y_val(1)],  'facecolor', [color.red 0.2], 'edgecolor', [color.red 0.2])
        hold off
        
        
        % Waveform
        try % see if the waveform file exists and works.
            load([this_S.label{1}(1:6) '-wv.mat'], 'mWV', 'xrange')
            
            subplot (3,4,1)
            hold on
            for ii  = 1:4
                plot (xrange(:,ii),mWV(:,ii),"LineWidth",3, 'color', c_ord(ii,:));
            end
%             legend("Ch 1", "Ch 2","Ch 3","Ch 4", 'fontsize', 8, 'location', 'north', 'orientation', 'horizontal' ); legend boxoff
            xlabel ("Time (ms)")
            axis off
            title ({strrep(fname, '_', ' ') ;  strrep(this_S.label{1}, '_', ' '); ['grade: ' this_S.label{1}(end-2)] ; Blocks{iB}})
            xlim([min(xrange,[], 'all') max(xrange,[], 'all')])
        catch
            subplot (3,4,1)
            text(0, .5, {'waveform file does not exist' ; 'or is corrupted'});
            axis off
        end
        
        % get the wave properties
        if exist('mWV', 'var')
            wave_prop = MS_get_wave_properties(this_S, [xrange(:,1) mWV],this_csc.tvec,  0);
            
        end
        
        % ISI
        subplot(3,4,2)
        histogram(diff(this_S.t{1})*1000,0:20:1000)
        xlabel('ISI (ms)')
        ylabel('spike count')
        xlim([0 1000]);
        y_val = ylim; 
        text(300, y_val(2)*.8, ['Fr: ' num2str(length(this_S.t{1}) / (this_csc.tvec(end) - this_csc.tvec(1))) 'Hz'])
        
        
        % Multiraster
        subplot (3,4,3:4)
        cfg_rast = [];
        cfg_rast.openNewFig = 0;
        yyaxis left
        MultiRaster(cfg_rast,this_S);
        ylim([0 1.5]); set(gca, 'ytick', []);
        ylabel([])
        hold on
        if iB == 2 || iB == 3
            yyaxis right
            plot(this_spd.tvec,this_spd.data)
            ylim([0 16])
            ylabel('speed (cm/s)');
%             legend({'spikes', 'speed'}, 'Orientation', 'horizontal', 'Location', 'northoutside'); legend boxoff;
        else
%             legend({'spikes'}, 'Orientation', 'horizontal', 'Location', 'northoutside'); legend boxoff;
        end
       
        
        
        
        %Restrict to periods of movement
        %         if iB == 2 || iB == 3
        %
        %             cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 5; % speed limit in cm/sec
        %             iv_fast = TSDtoIV(cfg,this_spd); % only keep intervals with speed above thresh
        %         end
        
        %     pos = restrict(pos,iv_fast);
        %     S = restrict(S,iv_fast);
        %
        % prepare S for plots
        if iB == 2 || iB == 3
            
            % interpolate the spikes to match the time vector
            move_idx =  (.5 < this_spd.data) & (this_spd.data < 6); 
            spk_x = interp1(this_pos.tvec(move_idx),this_pos.data(1,move_idx),this_S.t{1},'linear');
            spk_y = interp1(this_pos.tvec(move_idx),this_pos.data(2,move_idx),this_S.t{1},'linear');
            
            % plot
%             figure(101)
            % basic position plot.
            tvec_cord = winter(length(this_pos.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
            
            ax1 = subplot(3,4,5);
            scatter(this_pos.data(1,:), this_pos.data(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
            colormap(ax1, 'winter');
            cax = colorbar;
            cax.Position(1) = cax.Position(1) + 0.03;
            cax.Ticks = [0 1]; cax.TickLabels = {'start', 'end'};
            hold on
            
            % add some spikes
            %     S_idx = nearest_idx(S.t{1}, pos.tvec);
            
            plot(spk_x,spk_y, '.r', 'markersize', 6)
            axis off
            if iB == 3
                xlim([12 50]);
            end
            %% convert to heat map.
            
            
            % set up bins
            if iB == 2
            SET_xmin = 10; SET_ymin = 0; % set up bins
            SET_xmax = 82; SET_ymax = 60;
            SET_xBinSz = 2; SET_yBinSz =2;

            NaN_map = zeros(length(SET_xmin:SET_xBinSz:SET_xmax), length(SET_ymin:SET_yBinSz:SET_ymax));
            NaN_map(1:5,25:end) = NaN;
            NaN_map(30:end,25:end) = NaN;
            NaN_map(:,29:end) = NaN;
            
            nan_idx = isnan(NaN_map); 
            
            elseif iB == 3

            SET_xmin = 15; SET_ymin = 0; % set up bins
            SET_xmax = 45; SET_ymax = 30;
            SET_xBinSz = 1; SET_yBinSz =1;
            end
            
            
            x_edges = SET_xmin:SET_xBinSz:SET_xmax;
            y_edges = SET_ymin:SET_yBinSz:SET_ymax;
            
            % set up gaussian
            kernel = gausskernel([SET_xBinSz SET_yBinSz],SET_xBinSz/2); % 2d gaussian in bins
            
            
            % compute occupancy
            occ_hist = hist3(this_pos.data(1:2,move_idx)', 'edges', {x_edges y_edges});
            %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
            
            no_occ_idx = find(occ_hist < 1); % NaN out bins never visited
            occ_hist = conv2(occ_hist,kernel,'same');
            occ_hist(no_occ_idx) = NaN;
            if iB == 2; occ_hist(nan_idx) = NaN; end
            
            occ_hist = occ_hist .* mode(diff(this_pos.tvec)); % convert samples to seconds using video frame rate (30 Hz)
            
            subplot(3,4,6)
            pcolor(occ_hist'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'secs'; cb.Ticks = [0 cb.Ticks(end)];
            title('occupancy');
            
            % get the spike map
            spk_hist = hist3([spk_x, spk_y], 'edges', {x_edges y_edges});
            %     spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
            
            spk_hist = conv2(spk_hist,kernel,'same');
            spk_hist(no_occ_idx) = NaN;
            
            subplot(3,4,7)
            pcolor(spk_hist'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'nSpikes'; cb.Ticks = [0 cb.Ticks(end)];
            title('spikes');
            
            % rate map
            tc = spk_hist./occ_hist;
            
            subplot(3,4,8)
            pcolor(tc'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)];
            title('rate map');
            
            %% plot the heading by speed heatmap
%             cfg_speed = []; 
%             vector_mat = MS_speed_HD(cfg_speed,this_S, this_pos, this_spd);
%             
            
            %% split half
            Split.S{1} = restrict(S, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(1)+(block_idx.(Blocks{iB})(2) - block_idx.(Blocks{iB})(1))/2);
            Split.S{2} = restrict(S, block_idx.(Blocks{iB})(1)+(block_idx.(Blocks{iB})(2) - block_idx.(Blocks{iB})(1))/2,block_idx.(Blocks{iB})(2));
            
            Split.pos{1} = restrict(this_pos, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(1)+(block_idx.(Blocks{iB})(2) - block_idx.(Blocks{iB})(1))/2);
            Split.pos{2} = restrict(this_pos, block_idx.(Blocks{iB})(1)+(block_idx.(Blocks{iB})(2) - block_idx.(Blocks{iB})(1))/2,block_idx.(Blocks{iB})(2));
            
            Split.spd{1} = restrict(this_spd, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(1)+(block_idx.(Blocks{iB})(2) - block_idx.(Blocks{iB})(1))/2);
            Split.spd{2} = restrict(this_spd, block_idx.(Blocks{iB})(1)+(block_idx.(Blocks{iB})(2) - block_idx.(Blocks{iB})(1))/2,block_idx.(Blocks{iB})(2));

            for ii = 1:2 
            Split.move_idx{ii} =  (.5 < Split.spd{ii}.data) & (Split.spd{ii}.data < 6); 

            Split.spk_x{ii} = interp1(Split.pos{ii}.tvec(Split.move_idx{ii}),Split.pos{ii}.data(1,Split.move_idx{ii}),Split.S{ii}.t{1},'linear');
            Split.spk_y{ii} = interp1(Split.pos{ii}.tvec(Split.move_idx{ii}),Split.pos{ii}.data(2,Split.move_idx{ii}),Split.S{ii}.t{1},'linear');
            
            occ_hist = hist3(Split.pos{ii}.data(1:2,Split.move_idx{ii})', 'edges', {x_edges y_edges});
            %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
            
            no_occ_idx = find(occ_hist == 0); % NaN out bins never visited
            occ_hist = conv2(occ_hist,kernel,'same');
            occ_hist(no_occ_idx) = NaN;
            if iB == 2; occ_hist(nan_idx) = NaN; end

            occ_hist = occ_hist .* mode(diff(this_pos.tvec));
            
            spk_hist = hist3([Split.spk_x{ii}, Split.spk_y{ii}], 'edges', {x_edges y_edges});
            
            spk_hist = conv2(spk_hist,kernel,'same');
            spk_hist(no_occ_idx) = NaN;
            
            Split.tc{ii} = spk_hist./occ_hist;

                
            subplot(3,4,10+ii)
            pcolor(Split.tc{ii}'); shading flat; axis off; cb=colorbar; cb.Position(1) = cb.Position(1) + .03; cb.Label.String = 'rate (Hz)'; cb.Ticks = [0 cb.Ticks(end)];
            title(['TC half ' num2str(ii) ]);
            
            Split.occ_hist{ii} = occ_hist;
            Split.spk_hist{ii} = spk_hist;
            
            end
            idx = ~isnan(Split.tc{1}) & ~isnan(Split.tc{2});
            split_xcor = corr2(Split.tc{1}(idx), Split.tc{2}(idx));

            
        end
        
        % plot the autocorrelation if there is one
        if exist('wave_prop', 'var')
            subplot(3,4,10)
            bar(wave_prop.auto_corr.xbin*1000, wave_prop.auto_corr.ac)
            xlabel('Lag (ms)');
        end
        
        % get sta
%         cfg_sta = [];
%         [st_mat, t_win] = MS_get_sta(cfg_sta,this_csc,S);
%         subplot(3,4,12)
%         plot(t_win,nanmean(st_mat))
        
        
        % save the output
        SetFigure([], gcf);
        set(gcf, 'position', [81 -377  1760  880])
        
        saveas(gcf, [inter_dir filesep fname '_' S.label{1}(1:end-2) '_' Blocks{iB} '_check.png'])
        saveas(gcf, [inter_dir filesep fname '_' S.label{1}(1:end-2) '_' Blocks{iB} '_check.fig'])
        
        
        pause(1)
        close all
        
        %% export intermediate file
        This_Cell = [];
        This_Cell.Meta = Meta;
        This_Cell.ID = S.label{1};
        This_Cell.fname = fname;
        This_Cell.S = S;
        if exist('wave_prop', 'var')
            This_Cell.wave = wave_prop;
        else
            This_Cell.wave = [];
        end
        This_Cell.pos = pos;
        This_Cell.csc = csc;
        
        % block specifics. 
        
        This_Cell.(Blocks{iB}).psd.Px = Px;
        This_Cell.(Blocks{iB}).psd.Fx = Fx;
        if exist('wave_prop')
            This_Cell.(Blocks{iB}).wave = wave_prop;
        end
        if exist('st_mat', 'var')
            This_Cell.(Blocks{iB}).sta.mat = st_mat;
            This_Cell.(Blocks{iB}).sta.tvec = t_win;
        end
        if iB ==2 || iB ==3
            This_Cell.(Blocks{iB}).spatial.tc = tc;
            This_Cell.(Blocks{iB}).spatial.spk_hist = spk_hist;
            This_Cell.(Blocks{iB}).spatial.occ_hist = occ_hist;
            This_Cell.(Blocks{iB}).spatial.split_xcor = split_xcor;
        end

    end % end of the block
    
    mkdir([inter_dir filesep 'All_maze_cells']);
    save([inter_dir filesep 'All_maze_cells' filesep fname '_' S.label{1}(1:end-2) '.mat'],'This_Cell','-v7.3')
    
end % end of .t file list


