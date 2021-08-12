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
Meta = MS_Load_meta; 

% get position
cfg_pos = [];
cfg_pos.convFact = [8,8];
pos = MS_LoadPos(cfg_pos);
spd = getLinSpd([],pos);
% spd_gk = gausskernel(ceil(1/mode(diff(spd.tvec))),ceil(1/mode(diff(spd.tvec)))*.2);
% spd.data =  conv2(spd.data,gausswin(100),'same'); % smooth over 2 seconds
spd.data = smooth(spd.data, 2*ceil(1/mode(diff(spd.tvec))))'; 

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


block_idx = [];
if length(evt.t{start_idx}) > 1
    disp('splitting into blocks')
    Blocks = {'Pre_sleep', 'W_maze', 'OF', 'Post_sleep'};
    for iB = 1:length(Blocks)
    block_idx.(Blocks{iB}) = [evt.t{start_idx}(iB) evt.t{end_idx}(iB)]; 
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
    
    this_csc = [];
    this_csc = restrict(csc, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
    this_S = restrict(S, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
    this_pos = restrict(pos, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
    this_spd = restrict(spd, block_idx.(Blocks{iB})(1),block_idx.(Blocks{iB})(2));
    
    % generate psd
    cfg_psd.hann_win = 2^11;
    [Px, Fx] = pwelch(this_csc.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);

    figure(101)
    subplot(3,4,9:10)
    hold on
    plot(Fx, 10*log10(Px), 'color', [244,93,1]/255)
    xlim([0 150])
    set(gca,'XTick',0:10:150)
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
        load([this_S.label{1}(1:end-4) '-wv.mat'], 'mWV', 'xrange')
        
        subplot (3,4,1)
        hold on
        for ii  = 1:4
            plot (xrange(:,ii)+ii*2,mWV(:,ii),"LineWidth",3, 'color', c_ord(ii,:));
        end
        legend("Ch 1", "Ch 2","Ch 3","Ch 4", 'fontsize', 8); legend boxoff
        xlabel ("Time (ms)")
        axis off
        title ({strrep(fname, '_', ' ') ;  strrep(this_S.label{1}, '_', ' ')})
    catch
        subplot (3,4,1)
        text(0, .5, {'waveform file does not exist' ; 'or is corrupted'});
        axis off
    end
    
    % get the wave properties
    if exist('mv', 'var')
       wave_prop = MS_get_wave_properties(S, mv,pos.tvec,  0);  
        
    end
    
    % ISI
    figure(101)
    subplot(3,4,2)
    histogram(diff(this_S.t{1})*1000,0:20:1000)
    xlabel('ISI (ms)')
    ylabel('spike count')
    xlim([0 1000]); 
    
    % Multiraster
        figure(101)
    subplot (3,4,3:4)
    cfg_rast = [];
    cfg_rast.openNewFig = 0;
    yyaxis left
    MultiRaster(cfg_rast,this_S);
    ylim([0 1.5]); set(gca, 'ytick', []); 
    ylabel([])
    hold on
    yyaxis right
    plot(this_spd.tvec,this_spd.data)
    ylim([0 80])
    ylabel('speed (cm/s)');
    legend({'spikes', 'speed'}, 'Orientation', 'horizontal', 'Location', 'northoutside'); legend boxoff;
    
    
    
    
    %Restrict to periods of movement
    cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 5; % speed limit in cm/sec
    iv_fast = TSDtoIV(cfg,this_spd); % only keep intervals with speed above thresh
    
    %     pos = restrict(pos,iv_fast);
    %     S = restrict(S,iv_fast);
    %
    % prepare S for plots
    
    % interpolate the spikes to match the time vector
    spk_x = interp1(this_pos.tvec,this_pos.data(1,:),this_S.t{1},'linear');
    spk_y = interp1(this_pos.tvec,this_pos.data(2,:),this_S.t{1},'linear');
    
    % plot
    figure(101)
    % basic position plot.
    tvec_cord = winter(length(this_pos.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
    ax1 = subplot(3,4,5);
    scatter(this_pos.data(1,:), this_pos.data(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    
    colormap(ax1, 'winter');
    cax = colorbar;
    cax.Ticks = [0 1]; cax.TickLabels = {'start', 'end'};
    hold on
    
    % add some spikes
%     S_idx = nearest_idx(S.t{1}, pos.tvec);
    
    plot(spk_x,spk_y, '.r', 'markersize', 12)
    axis off
    %% convert to heat map.
    
    
    % set up bins
    SET_xmin = 0; SET_ymin = 10; % set up bins
    SET_xmax = 60; SET_ymax = 70;
    SET_xBinSz = 2; SET_yBinSz = 2;
    
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    % set up gaussian
    kernel = gausskernel([SET_xBinSz/2 SET_yBinSz/2],SET_xBinSz/2); % 2d gaussian in bins
    
    
    % compute occupancy
    occ_hist = hist3(this_pos.data(1:2,:)', 'edges', {y_edges x_edges});
    %     occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    
    occ_hist = conv2(occ_hist,kernel,'same');
    
    no_occ_idx = find(occ_hist == 0); % NaN out bins never visited
    occ_hist(no_occ_idx) = NaN;
    
    occ_hist = occ_hist .* (1/30); % convert samples to seconds using video frame rate (30 Hz)
    
    subplot(3,4,6)
    pcolor(occ_hist'); shading flat; axis off; colorbar
    title('occupancy');
    
    % get the spike map
    spk_hist = hist3([spk_x, spk_y], 'edges', {y_edges x_edges});
    %     spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    
    spk_hist = conv2(spk_hist,kernel,'same');
    
    
    spk_hist(no_occ_idx) = NaN;
    
    subplot(3,4,7)
    pcolor(spk_hist'); shading flat; axis off; colorbar
    title('spikes');
    
    % rate map
    tc = spk_hist./occ_hist;
    
    subplot(3,4,8)
    pcolor(tc'); shading flat; axis off; colorbar
    title('rate map');
    
    
    %% plot the autocorrelation if there is one
    if exist('wave_prop', 'var')
        subplot(3,4,11)
        bar(wave_prop.auto_corr.xbin*1000, wave_prop.auto_corr.ac)
        xlabel('Lag (ms)');
    end
    
    %% get sta
    cfg_sta = [];
    [st_mat, t_win] = MS_get_sta(cfg_sta,csc,S);
    subplot(3,4,12)
    plot(t_win,nanmean(st_mat))
    
    
    %% save the output
    SetFigure([], gcf);
    set(gcf, 'position', [81 -377  1760  880])
    
    saveas(gcf, [inter_dir filesep fname '_check.png'])
    saveas(gcf, [inter_dir filesep fname '_check.fig'])

    
    pause(1)
    %close all
    
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
    This_Cell.psd.Px = Px;
    This_Cell.psd.Fx = Fx;
    This_Cell.sta.mat = st_mat;
    This_Cell.sta.tvec = t_win;
    This_Cell.spatial.tc = tc;
    This_Cell.spatial.spk_hist = spk_hist;
    This_Cell.spatial.occ_hist = occ_hist;

    mkdir([inter_dir filesep 'All_cells']); 
    save([inter_dir filesep 'All_cells' filesep fname '_' S.label{1}(1:end-2) '.mat'],'This_Cell','-v7.3')
end % end of .t file list