function out = MS_get_place_field(cfg_in, S, pos)
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


cfg_def = [];




cfg = ProcessConfig(cfg_def, cfg_in); 


%% set up bins
SET_xmin = floor(min(pos.data(2,:))/100)*100 ; SET_ymin = floor(min(pos.data(1,:))/100)*100; % set up bins
SET_xmax = ceil(max(pos.data(2,:))/100)*100;   SET_ymax = ceil(max(pos.data(1,:))/100)*100;

% SET_xmin = 100; 
% SET_xmax = 400;
% SET_ymin = 200;
% SET_ymax = 700;
SET_xBinSz =30; SET_yBinSz = 30;


 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([2 2],2); % 2d gaussian in bins

% compute occupancy encode
occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');
 
no_occ_idx = find(occ_hist < 2); % NaN out bins never visited
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
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx) = NaN;
    tc = spk_hist./occ_hist;

    subplot(m, n, s2_idx(ip))
    pcolor(tc'); shading flat; axis off; %colorbar('Location', 'northoutside')
   
end