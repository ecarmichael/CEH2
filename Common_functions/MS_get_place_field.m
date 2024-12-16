function out = MS_get_place_field(cfg_in, S, pos, spd)
%% MS_get_place_field: computes the rate map for a given set of cells.
%
%
%
%    Inputs:
%    - cfg_in: struct]    configurations
%
%    - S [struct]   TS struct containing spike times
%
%    - pos [struct]    TSD for the position.
%
%    - spd [struct]  tsd for speed. Can be used to limit to periods of
%    movement. if empty will not limit.
%
%    Outputs:
%    - out: [1 x NCell]   spatial information and rate map for each cell.
%
%
%
%
% EC 2023-07-10   initial version
%
%
%
%% initialize

if nargin < 4
    spd = [];
end


cfg_def = [];
cfg_def.gau_sd = 1;
cfg_def.gau_win =2;
cfg_def.xBinSz = 2.5;
cfg_def.yBinSz = 2.5;
cfg_def.min_occ = .25;
cfg_def.plot = 1;
cfg_def.spd_thresh = 2.5; % in cm/s;

cfg = ProcessConfig(cfg_def, cfg_in);

% gaussian kernal
kernel = gausskernel([cfg.gau_win cfg.gau_win],cfg.gau_sd); % 2d gaussian in bins

%% restrict to movement if speed is given.

if ~isempty(spd) || sum(contains(pos.label, 'Speed'))>0
    
    
    
    spd = pos;
    this_idx = find(contains(pos.label, 'Speed'));
    spd.label(~this_idx) = [];
    spd.data = spd.data(this_idx,:);
    
    if ~isstruct(spd) && length(spd) == length(pos.tvec);
        
        spd = pos;
        spd.data = spd;
        
        
    end
    
    cfg_spd = [];
    cfg_spd.method = 'raw';
    cfg_spd.operation = '>';
    cfg_spd.threshold = cfg.spd_thresh;
    iv_move = TSDtoIV(cfg_spd,spd); % only keep intervals with speed above thresh
    
    
    S = restrict(S, iv_move);
    pos = restrict(pos, iv_move);
    
end


%% set up bins
cfg.xmin = floor(min(pos.data(1,:))/10)*10 ; cfg.ymin = floor(min(pos.data(2,:))/10)*10; % set up bins
cfg.xmax = ceil(max(pos.data(1,:))/10)*10;   cfg.ymax = ceil(max(pos.data(2,:))/10)*10;

% cfg.xmin = 100;
% cfg.xmax = 400;
% cfg.ymin = 200;
% cfg.ymax = 700;
% cfg.xBinSz =30; cfg.yBinSz = 30;



x_edges = cfg.xmin:cfg.xBinSz:cfg.xmax;
y_edges = cfg.ymin:cfg.yBinSz:cfg.ymax;


% compute occupancy encode
occ_hist = histcn(pos.data(1:2,:)',x_edges,y_edges); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');


occ_hist = occ_hist .* (1/30); % convert samples to seconds using video frame rate (30 Hz)
no_occ_idx = find(occ_hist < cfg.min_occ); % NaN out bins never visited
occ_hist(no_occ_idx) = NaN;


%% compute the tunig curves
for ii = 1:length(S.t)
    
    if isempty(S.t{ii})
        tc{ii} = NaN(size(occ_hist)); 
    else
    
    
    spk_x{ii} = interp1(pos.tvec,pos.data(1,:),S.t{ii},'linear');
    spk_y{ii} = interp1(pos.tvec,pos.data(2,:),S.t{ii},'linear');
    
    
    % compute the rate map
    spk_hist = histcn([spk_x{ii}, spk_y{ii}],x_edges,y_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx) = NaN;
    tc{ii} = spk_hist./occ_hist;
    end
    
end
%% plot some simple maps with spikes

if cfg.plot
    m = 4;
    n = 4;
    s1_idx = 1:2:n*m;
    s2_idx = 2:2:n*m;
    figure(102)
    clf
    ip = 0;
    for ii = 1:length(S.t)
        
        if isempty(S.t{ii})
            continue
        end
        
        if ip >= (n*m)/2
            figure(get(gcf, 'Number')+1)
            ip = 1;
        else
            ip = ip+1;
        end
        
        disp(ip)
        subplot(m, n, s1_idx(ip))
        plot(pos.data(1,:), pos.data(2,:), '.k');
        hold on
        plot(spk_x{ii}, spk_y{ii}, '.r');axis off;
        
        
        
        subplot(m, n, s2_idx(ip))
        pcolor(tc{ii}'); shading flat; axis off; %colorbar('Location', 'northoutside')
        
        title([strrep(S.label{ii}(1:end-2), '_', '-' ) ' (peak: ' num2str(max(tc{ii},[], 'all'),2) 'Hz)'])
        
    end
end
%%

out = [];
out.label = S.label;
out.cfg = cfg;
out.tc = tc;
out.occ = occ_hist;
