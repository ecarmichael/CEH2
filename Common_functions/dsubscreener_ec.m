function dsubscreener_ec(data_dir, inter_dir)
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

%% get all .t files

file_list = FindFiles('*.t');





for iC = 1:length(file_list)
    
    cfg = [];
    cfg.fc = file_list(iC);
    cfg.getTTnumbers = 0;
    
    S = LoadSpikes(cfg);
    
    % generate psd
    cfg_psd.hann_win = 2^11;
    [Px, Fx] = pwelch(csc.data, hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , csc.cfg.hdr{1}.SamplingFrequency);

    figure(101)
    subplot(3,4,9:10)
    hold on
    plot(Fx, 10*log10(Px))
    xlim([0 150])
    set(gca,'XTick',0:10:150)
    xlabel("Frequency")
    hold off


    % Waveform
    try % see if the waveform file exists and works. 
        mv= csvread([S.label{1}(1:3) '_AvgWaveforms.csv']);
        
        subplot (3,4,1)
        hold on
        for ii  = 1:4
            plot (mv(1,:)+ii*2,mv(ii+1,:),"LineWidth",3, 'color', c_ord(ii,:));
        end
        legend("Ch 1", "Ch 2","Ch 3","Ch 4", 'fontsize', 8); legend boxoff
        xlabel ("Time (ms)")
        axis off
        title ({strrep(fname, '_', ' ') ;  strrep(S.label{1}, '_', ' ')})
    catch 
                subplot (3,4,1)
    text(0, .5, {'waveform file does not exist' ; 'or is corrupeted'});
    axis off 
    end
    
    % ISI
    subplot(3,4,2)
    histogram(diff(S.t{1})*1000,0:20:1000)
    xlabel('ISI (ms)')
    ylabel('spike count')
    xlim([0 1000]); 
    
    % Multiraster
    
    subplot (3,4,3:4)
    cfg_rast = [];
    cfg_rast.openNewFig = 0;
    yyaxis left
    MultiRaster(cfg_rast,S);
    ylim([0 1.5]); set(gca, 'ytick', []); 
    ylabel([])
    hold on
    yyaxis right
    plot(spd.tvec,spd.data)
    ylim([0 80])
    ylabel('speed (cm/s)');
    legend({'spikes', 'speed'}, 'Orientation', 'horizontal', 'Location', 'northoutside'); legend boxoff;
    
    
    
    
    
    %Restrict to periods of movement
    cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 5; % speed limit in cm/sec
    iv_fast = TSDtoIV(cfg,spd); % only keep intervals with speed above thresh
    
    %     pos = restrict(pos,iv_fast);
    %     S = restrict(S,iv_fast);
    %
    % prepare S for plots
    
    % interpolate the spikes to match the time vector
    spk_x = interp1(pos.tvec,pos.data(1,:),S.t{1},'linear');
    spk_y = interp1(pos.tvec,pos.data(2,:),S.t{1},'linear');
    
    % plot
    figure(101)
    % basic position plot.
    tvec_cord = winter(length(pos.data(1,:)));% repmat(0.2, length(pos.data(1,:)),1)];
    ax1 = subplot(3,4,5);
    scatter(pos.data(1,:), pos.data(2,:), 55, tvec_cord, '.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    
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
    occ_hist = hist3(pos.data(1:2,:)', 'edges', {y_edges x_edges});
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
    
    
    %% save the output
    SetFigure([], gcf);
    set(gcf, 'position', [81 -377  1760  880])
    
    saveas(gcf, [inter_dir filesep fname '_check.png'])
    saveas(gcf, [inter_dir filesep fname '_check.fig'])

    
    pause(1)
    close all
end % end of .t file list