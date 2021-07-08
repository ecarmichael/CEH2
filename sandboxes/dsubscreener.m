function dsubscreener(data_dir)
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


%% get all .t files

file_list = FindFiles('*.t');

   % get position
    cfg.getTTnumbers = 0;
    S = LoadSpikes(cfg);


    pos = MS_LoadPos([]);

   % get events
    spk_x = interp1(pos.tvec,pos.data(1,:),S.t{1},'linear');
    spk_y = interp1(pos.tvec,pos.data(2,:),S.t{1},'linear');
   


for iC = 1:length(file_list)
        
    cfg = [];
    cfg.fc = {file_list{iC}};
    cfg.getTTnumbers = 0;
    
    
    S = LoadSpikes(cfg);
    
    
    % prepare S for plots
    
    % interpolate the spikes to match the time vector
    spk_x = interp1(pos.tvec,pos.data(1,:),S.t{1},'linear');
    spk_y = interp1(pos.tvec,pos.data(2,:),S.t{1},'linear');
    
    % plot
    figure(101)
    
    subplot(3,3,5)
    plot(pos.data(1,:), pos.data(2,:), '.', 'color', [0.8 0.8 0.8]);
    hold on
    
    
    S_idx = nearest_idx(S.t{1}, pos.tvec);
    
    
    plot(spk_x,spk_y, '.r')
    axis off
    %% convert to heat map.
    
    % set up bins
    SET_xmin = 80; SET_ymin = 0; % set up bins
    SET_xmax = 660; SET_ymax = 520;
    SET_xBinSz = 20; SET_yBinSz = 20;
    
    
    x_edges = SET_xmin:SET_xBinSz:SET_xmax;
    y_edges = SET_ymin:SET_yBinSz:SET_ymax;
    
    
    % compute occupancy
    occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
    
    no_occ_idx = find(occ_hist == 0); % NaN out bins never visited
    occ_hist(no_occ_idx) = NaN;
    
    occ_hist = occ_hist .* (1/30); % convert samples to seconds using video frame rate (30 Hz)
    
    subplot(3,3,6)
    pcolor(occ_hist); shading flat; axis off; colorbar
    title('occupancy');
    
    % get the spike map
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    
    spk_hist(no_occ_idx) = NaN;
    
    subplot(3,3,8)
    pcolor(spk_hist); shading flat; axis off; colorbar
    title('spikes');
    
    % rate map
    tc = spk_hist./occ_hist;
    
    subplot(3,3,9)
    pcolor(tc); shading flat; axis off; colorbar
    title('rate map');
    
    % Waveform
    
    mv= [csvread("TT3_AvgWaveforms.csv")];
    
    subplot (3,3,4)
    plot (mv(1,:),mv(2:end,:),"LineWidth",3)
    legend("Channel 1", "Channel 2","Channel 3","Channel 4")
    xlabel ("Time (ms)")
    title ("M23 Tetrode 3 July 2nd")
    
    
    %% ISI
    subplot(3,3,7)
    histogram(diff(S.t{1})*1000,0:20:1000)
    
    %% Label
    subplot (3,3,1)
    plot(nan(1,10))
    text(0,.8,["Cell:" strrep(S.label{1},"_"," ")],"fontsize",20)
    
    %% Multiraster
    cfg.openNewFig==0
    subplot (3,3,2:3)
    MultiRaster([],S)
    
end % end of .t file list