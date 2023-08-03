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
cfg_def.min_occ = 7.5; 

cfg = ProcessConfig(cfg_def, cfg_in); 

% gaussian kernal
kernel = gausskernel([cfg.gau_win cfg.gau_win],cfg.gau_sd); % 2d gaussian in bins

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

no_occ_idx = find(occ_hist < cfg.min_occ); % NaN out bins never visited
occ_hist(no_occ_idx) = NaN;
 

 
occ_hist = occ_hist .* (1/30); % convert samples to seconds using video frame rate (30 Hz)

%% plot some simple maps with spikes
m = 4;
n = 4;
s1_idx = 1:2:n*m; 
s2_idx = 2:2:n*m; 
figure(102)
clf
ip = 0; 
tc = [];
for ii = 1:length(S.t)
    spk_x = interp1(pos.tvec,pos.data(1,:),S.t{ii},'linear');
    spk_y = interp1(pos.tvec,pos.data(2,:),S.t{ii},'linear');
    if ip >= (n*m)/2
        figure(102+ii)
        ip = 1;
    else
        ip = ip+1;
    end
    
    disp(ip)
    subplot(m, n, s1_idx(ip))
    plot(pos.data(1,:), pos.data(2,:), '.k');
    title(S.label{ii})
    hold on
    plot(spk_x, spk_y, '.r');axis off;
    
    % compute the rate map
    spk_hist = histcn([spk_x, spk_y],x_edges,y_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx) = NaN;
    tc{ii} = spk_hist./occ_hist;

    subplot(m, n, s2_idx(ip))
    pcolor(tc{ii}'); shading flat; axis off; %colorbar('Location', 'northoutside')
   
end
%% 

out = [];
out.cfg = cfg; 
out.tc = tc; 
out.occ = occ_hist;
