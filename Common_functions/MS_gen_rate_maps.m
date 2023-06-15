function MS_gen_rate_maps
%% MS_gen_rate_maps: loads all .t files and position in a dir and generates basic rate maps for screening. Hacky, but helpful. 







%% load the pos and S
evt = LoadEvents([]);
start_t = evt.t{1}(1); 

cfg = [];
cfg.fc = {'CSC1.ncs'};
csc = LoadCSC(cfg); 
    fprintf('CSC             %2.2d    |    %2.2d \n', csc.tvec(1), csc.tvec(end) - csc.tvec(1))

pos = LoadPos([]);
% pos.tvec = pos.tvec - start_t; 
    fprintf('Pos             %2.2d    |    %2.2d \n', pos.tvec(1), pos.tvec(end) - pos.tvec(1))

cfg = [];
cfg.getTTnumbers = 0; 
cfg.uint = '64'; 
% cfg.tsflag = 'ts'; 
cfg.verbose = 0; 
S = LoadSpikes(cfg);

for ii = 1:length(S.label)
    fprintf('S %s   %2.2d    |    %2.2d \n',S.label{ii}, (S.t{ii}(1)), S.t{ii}(end) - S.t{ii}(1))
end
    





%%  get the occupancy



% set up bins
SET_xmin = floor(min(pos.data(2,:))/100)*100 ; SET_ymin = floor(min(pos.data(1,:))/100)*100; % set up bins
SET_xmax = ceil(max(pos.data(2,:))/100)*100; SET_ymax = ceil(max(pos.data(1,:))/100)*100;

% SET_xmin = 100; 
% SET_xmax = 400;
% SET_ymin = 200;
% SET_ymax = 700;
SET_xBinSz =15; SET_yBinSz = 15;


 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([1 1],1); % 2d gaussian in bins

% compute occupancy encode
occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');
 
no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
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