%% EMD border score sandbox


cd('/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/ck2cre-1359hd/2021_01_30/14_18_06')
load('behav_DLC.mat', 'behav')
load('ms.mat', 'ms')

%% get ms binary

ms = MS_msExtractBinary_detrendTraces(ms,2);

%% exclude non-movement periods

% check if the behav needs to be interpolated.
if behav.time(end) ~= ms.time(end) || length(behav.time) ~= length(ms.time)
    fprintf('<strong> %s </strong>: behaviour and Ca are not the same length or end time.  attempting alignment \n', mfilename);
    behav_aligned = MS_align_data(behav, ms);
else
    behav_aligned = behav;
end

%smooth speed
if exist('smooth', 'builtin')
    behav_aligned.speed = smooth(behav_aligned.speed, 3*mode(diff(ms.time)));
else
    
    behav_aligned.speed= conv2(behav_aligned.speed, gausswin(3*mode(diff(ms.time))), 'same');
end

movement_idx = behav_aligned.speed >2.5; % get times when the animal was moving.
accel_movement_idx = behav_aligned.speed(1:end-1) >2.5; % same as above but correct for diff used in acceleration calculation.

% get the acceleration
behav_aligned.accel = diff(behav_aligned.speed);


% compute event histogram
evts_idx = ms.Binary(:,1) == 1;


%% build some bins

cfg.p_thres = 0.05; % value for pvalue cut off;
cfg.stability_thres = 0.5; % from van der Veldt 2020
cfg.nShuff = 600;
cfg.p_bin_size = 3 ; % in cm
cfg.split_gaus_sd = 3; % sd for gaussian smoothing of place tuning for split session xcorr.
cfg.gauss_sd = 1; 
cfg.gauss_f_size = 3;


X_bins = 0:cfg.p_bin_size:ceil(max(behav_aligned.position(:,1)));
X_bin_centers = X_bins +  cfg.p_bin_size/2;
X_bin_centers = X_bin_centers(1:end-1);
% same for Y bins
Y_bins = 0:cfg.p_bin_size:ceil(max(behav_aligned.position(:,2)));
Y_bin_centers = Y_bins +  cfg.p_bin_size/2;
Y_bin_centers = Y_bin_centers(1:end-1);

%% generate some templates
z_temp = zeros(length(X_bins), length(Y_bins));
templates = [];

templates.B = z_temp;
templates.B(1:end,1) = 1; templates.B(1:end,end) = 1;
templates.B(1,1:end) = 1; templates.B(end,1:end) = 1;

%north
templates.N = z_temp;
templates.N(1,1:end) = 1;

%East
templates.E = z_temp;
templates.E(1:end,end) = 1;

%south
templates.S = z_temp;
templates.S(end,1:end) = 1;

%East
templates.W = z_temp;
templates.W(1:end,1) = 1;

temps = fieldnames(templates);


% smooth the templates and plot them 

figure(401)
set(gcf, 'position', [400  725 1100 160])
for ii = 1:length(temps)
    
%         2D guass filter.  No distoritions like kernal conv on E and S walls. 
     templates.(temps{ii}) = imgaussfilt( templates.(temps{ii}),cfg.gauss_sd, 'FilterSize',cfg.gauss_f_size);


    subplot(1,length(temps),ii);
       imagesc(templates.(temps{ii}))
    xlabel(temps{ii})
end

% % square the subplots.  Ugly but at least it's sqaure
% axesHandles = findobj(get(401,'Children'), 'flat','Type','axes');
% % Set the axis property to square
% axis(axesHandles,'square')
%% plot an example
iC  = 1;
if ishandle(101)
    close(101)
end
figure(101)
subplot(4,3,1:4)
hold on
plot(ms.time, ms.RawTraces(:,iC));
plot(ms.time, ms.deconvolvedSig(:,iC));
plot(ms.time, ms.denoisedCa(:,iC));
plot(ms.time, ms.Binary(:,iC));
legend('Raw', 'Decon', 'Denoised', 'Binary', 'Orientation', 'horizontal','Location','North');
title(['Cell:' num2str(iC)])


% subplot(4,3,[5 6 9 10])


%% generate a post prob map as a proxy for a rate map.




% get place information
[Place_MI, Place_posterior, Place_occupancy, ~, Place_tuning_curve] = MS_get_spatial_information_2D(ms.Binary(movement_idx,iC),behav_aligned.position(movement_idx,:), X_bins, Y_bins );

% get shuffle data
Place_shuff_tuning_curve = MS_split_shuff(ms.Binary(:,iC), behav_aligned.position,movement_idx, cfg.nShuff, X_bins, Y_bins);

% get stats
Place_stats= MS_boot_shuff(ms.Binary(:,iC), behav_aligned.position,movement_idx, cfg.nShuff, X_bins, Y_bins);

%get sig tuning
pval = sum(Place_shuff_tuning_curve > Place_tuning_curve,3)/cfg.nShuff;
Place_Sig_map = Place_tuning_curve;
Place_Sig_map(pval > cfg.p_thres) = 0;

%% generate a


%% check for sig using mean vector length vs shuffle
nShuff = 10000;
m_idx = find(movement_idx == 1);
e_idx = find(evts_idx == 1);

for iShuff = nShuff:-1:1
    
    shuff_idx = m_idx(randperm(length(m_idx),length(e_idx)));
    
    shuf_r(iShuff) = circ_r(behav_aligned.HD(shuff_idx));
    
end

p95_shuff_r = prctile(shuf_r, 95);

evt_r = circ_r(behav_aligned.HD(movement_idx & evts_idx));
occ_r = circ_r(behav_aligned.HD(movement_idx));
%
% hold on
% h= polarhistogram(behav_aligned.HD(movement_idx), 15);
% h.DisplayStyle = 'stairs';
% h.EdgeColor = [.7 .7 .7];
% h.EdgeAlpha = .4;
%
% h= polarhistogram(behav_aligned.HD(movement_idx & evts_idx), 15);
% h.DisplayStyle = 'stairs';
% h.EdgeColor = 'b';


% try with histograms and polar plots
figure(102)
% subplot(1,2,1)
% hold on
%
% [Occ, bins] = histcounts(behav_aligned.HD(movement_idx), -180:15:180);
% bin_c = bins(1:end-1)+15/2;
% bin_r = deg2rad(bin_c);
%
% Occ_norm = Occ/max(Occ);
% bar(bin_r, Occ_norm, 'facecolor', [.7 .7 .7], 'facealpha', .4)
%
% [Act, ~] = histcounts(behav_aligned.HD(movement_idx & evts_idx), -180:15:180);
% Act_norm = Act/max(Act);
%
%
% bar(bin_r, Act_norm, 'r')
% set(gca, 'xtick', [-pi:pi/2:pi], 'xticklabel', {'-pi', '-1/2pi', '0', '1/2pi', 'pi'})
%
% subplot(1,2,2)
hold on

polarhistogram(behav_aligned.HD(movement_idx), 15, 'Normalization','probability', 'DisplayStyle','stairs', 'edgecolor', [0.7 0.7 0.7], 'linewidth', 4)
% text(gca,num2str(occ_r))
polarhistogram(behav_aligned.HD(movement_idx & evts_idx), 15, 'Normalization','probability', 'DisplayStyle','stairs', 'edgecolor', 'b', 'linewidth', 4)




%% simple plot for a few cells
figure(101)

% subplot(3,2,1)
h = polarplot(behav_aligned.HD(movement_idx),ms.Binary(movement_idx));
% h.DisplayStyle = 'stairs';


%% gen some fake cell


HD_cell = zeros(1, 10);

% b_cell(1,2:8) = [2:2:6, 8,9, 4:-2:2];
HD_cell(1,2:8) = [2:2:6, 8,9, 4:-2:2]*4;
HD_cell = HD_cell./max(HD_cell,[],'all'); % normalize activity.
% b_cell = imgaussfilt(b_cell, 1);

figure(20)
subplot(2,2,1)
imagesc(HD_cell)


kernel_size = [size(HD_cell,1) size(HD_cell,2)];
occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
kernel = pdf('Normal', Rgrid, 0, occupancy_std);
kernel = kernel./sum(sum(kernel));

b_cell_s = conv2(HD_cell, kernel, 'same');

b_cell_s(b_cell_s < max(b_cell_s, [], 'all')*.3) = 0;

subplot(2,2,3)
imagesc(b_cell_s)

% generate a center cell for comparison
c_cell = zeros(10, 10);

c_cell(5:7,4:6) = [4,6,4; 6, 8, 6; 4, 6, 4];
c_cell = c_cell./max(c_cell,[],'all'); % normalize activity.


kernel_size = [size(c_cell,1) size(c_cell,2)];
occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
kernel = pdf('Normal', Rgrid, 0, occupancy_std);
kernel = kernel./sum(sum(kernel));

c_cell_s = conv2(c_cell, kernel, 'same');

c_cell_s(c_cell_s < max(c_cell_s, [], 'all')*.3) = 0;


% generate a corner edge cell for comparison
e_cell = zeros(10, 10);

e_cell(1:3,1:3) = [4,6,4; 6, 8, 6; 4, 6, 4];
e_cell = e_cell./max(e_cell,[],'all'); % normalize activity.


kernel_size = [size(e_cell,1) size(e_cell,2)];
occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
kernel = pdf('Normal', Rgrid, 0, occupancy_std);
kernel = kernel./sum(sum(kernel));

e_cell_s = conv2(e_cell, kernel, 'same');

e_cell_s(e_cell_s < max(e_cell_s, [], 'all')*.3) = 0;


%% compute the border score based on solstad 2008

% for testing
wall_vals = reshape(1:numel(b_cell_s), size(b_cell_s))';


wall_id = {'N', 'E', 'S', 'W'};
wall_score(1,:) =  sum(b_cell_s(1, 1:size(b_cell_s,2)) ~= 0)/size(b_cell_s,2);
wall_score(2,:) =  sum(b_cell_s(1:size(b_cell_s,2),size(b_cell_s, 1)) ~= 0)/size(b_cell_s,1);
wall_score(3,:) =  sum(b_cell_s(size(b_cell_s,2),1:size(b_cell_s,1)) ~= 0)/size(b_cell_s,1);
wall_score(4,:)=  sum(b_cell_s(1:size(b_cell_s,2),size(b_cell_s, 1)) ~= 0)/size(b_cell_s,2);

% get the all indicies
wall_idx(1,:) = wall_vals(1, 1:size(b_cell_s,2));
wall_idx(2,:) = wall_vals(1:size(b_cell_s,2),size(b_cell_s, 1));
wall_idx(3,:) = wall_vals(size(b_cell_s,2),1:size(b_cell_s,1));
wall_idx(4,:) = wall_vals(1:size(b_cell_s,2),size(b_cell_s, 1));

[cM, idx] = max(wall_score);

% get the distance from the wall
act_bins = find((HD_cell ~=0)');
[act_i, act_j] = find((HD_cell ~=0)');

range_mat = reshape(1:(size(HD_cell,1) * size(HD_cell,2)), size(HD_cell))';

peak = max(HD_cell,[], 'all');
[peak_idx_i, peak_idx_j] = find(HD_cell == peak);

dist = [];
for iB = length(act_bins):-1:1 % look across active pixels
    for iW = size(wall_idx,1):-1:1
        if contains(wall_id(iW), {'N' 'S'})
            [ii, jj] = find(wall_vals == wall_idx(iW,:));
        else
            [ii, jj] = find(wall_vals == wall_idx(iW,:)');
        end
        [ii_pt, jj_pt] = find(range_mat == act_bins(iB));
        dist(:,iB, iW) = sqrt(((ii - ii_pt).^2)+(jj - jj_pt).^2);
    end
end

% get the distance to the closest wall

dist_m = min(dist,[], 1);

[~, best_wall] = min(sum(dist_m,2));
fprintf('Best wall is <strong>%s</strong>\n', wall_id{best_wall})

% norm to largest possible distance
dM = mean(squeeze(dist_m(:,:,best_wall)./(min(size(b_cell_s))/2))); % normalize to max wall length and get mean.


b_score = (cM-dM)/(cM+dM);
fprintf('Border score = %d</strong>\n', b_score)


%% check plot

figure(20)
subplot(2,3,1)
imagesc(HD_cell)
b_score = MS_border_score(HD_cell);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,2)
imagesc(c_cell)
b_score = MS_border_score(c_cell);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,3)
imagesc(e_cell)
b_score = MS_border_score(e_cell);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,4)
imagesc(b_cell_s)
b_score = MS_border_score(b_cell_s);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,5)
imagesc(c_cell_s)
b_score = MS_border_score(c_cell_s);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,6)
imagesc(e_cell_s)
b_score = MS_border_score(e_cell_s);
text(0, 0, ['Border score = ' num2str(b_score,2)])

%% Get borderscores for all cells

for iC = 1:length(All_cells.place.MI)
    b_score = MS_border_score(abs(All_cells.place.Sig_map(:,:,iC)), 1, 1);
    All_cells.border.score(iC) = b_score;
end


figure(303)
bins = -1:.1:1;
h = histogram(All_cells.border.score, bins) ;
h.FaceColor = PARAMS.red;  h.FaceAlpha = .9;
h.EdgeColor = 'k';  h.EdgeAlpha = 1;
h.LineWidth = 1;
xlabel('border score')
ylabel('count')
title('border score for ck2cre-1359 2021-01-30')
vline(0.5, 'r', 'thresh')

for iC = 1:round(length(All_cells.border.score)/3)
    if All_cells.border.score(iC) >= 0.4  && sum(All_cells.place.Sig_map(:,:,iC) ~= 0, 'all') > 3
        [~,~,~,this_mat] = MS_border_score(All_cells.place.Sig_map(:,:,iC), 1, 1);
        
        figure(1000+ iC)
        imagesc(0:3:30, 0:3:30, this_mat)
        axis xy
        text(0,33, ['Cell ID: ' num2str(iC) '  Border Score: ' num2str(All_cells.border.score(iC), 2)])
    end
end

