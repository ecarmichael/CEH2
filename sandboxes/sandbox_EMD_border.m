%% EMD border score sandbox
addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2')); 
addpath('/home/ecarmichael/Documents/GitHub/OASIS_matlab');
oasis_setup;
% addpath(genpath('/home/ecarmichael/Documents/GitHub/EMDeBruijn/FastEMD')); %Add the EMD repo
% cd('/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/ck2cre-1359hd/2021_01_30/14_18_06')
% load('behav_DLC.mat', 'behav')
% data_dir = '/mnt/Data/Williams_Lab/II_classification/msbehavplace/ck2cre1/8-24-20/ck2cre1_sqOF_H14_M53_S37';
% data_dir = '/mnt/Data/Williams_Lab/II_classification/msbehavplace/ck2cre1/8-26-20/ck2cre1_sqOF-H14_M20_S21';
data_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/ck2cre-1359hd/2021_01_30/14_18_06'; 
parts = strsplit(data_dir, filesep); 

% save_dir = ['/mnt/Data/Williams_Lab/II_classification/Inter' filesep parts{end-2} filesep datestr(parts{end-1}, 'yyyy-mm-dd') filesep parts{end}(9:10) filesep 'Border'];
% mkdir(save_dir);
cd(data_dir) 

if exist('behav_DLC.mat','file')
    load('behav_DLC.mat');
else
load('behav.mat');
end
load('ms.mat', 'ms')
%% align behav

if behav.time(end) ~= ms.time(end) || length(behav.time) ~= length(ms.time)
    fprintf('<strong> %s </strong>: behaviour and Ca are not the same length or end time.  attempting alignment \n', mfilename);
    behav_aligned = MS_align_data(behav, ms);
else
    behav_aligned = behav;
end

split_idx = ceil(length(behav_aligned.time)/2); 
move_idx = behav_aligned.speed >2; % get times when the animal was moving.
move_idx_1 = behav_aligned.speed(1:split_idx)>2;
move_idx_2= behav_aligned.speed(split_idx:end)>2;

%% decon

ms = MS_msExtractBinary_detrendTraces(ms);

for iC = size(ms.RawTraces,2):-1:1
    ms.detrendRaw(:,iC) = detrend(ms.RawTraces(:,iC),2);
    [denoise(:,iC),deconv(:,iC)] = deconvolveCa(ms.detrendRaw(:,iC), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
    
end

%% generate a rate map
bin_size = 3;
plot_flag = 0;
smooth_sd = 2; %sd for smoothing

X_bins = 0:bin_size:ceil(max(behav_aligned.position(:,1)));
Y_bins = 0:bin_size:ceil(max(behav_aligned.position(:,2)));


for iC = size(ms.RawTraces,2):-1:1
    [rate_map(:,:,iC), occ_map(:,:,iC)] = MS_decon_rate_map(deconv(move_idx,iC), ms.time(move_idx), behav_aligned.position(move_idx,:), X_bins,Y_bins, plot_flag, smooth_sd);
    
    if plot_flag
        subplot(2,2,1)
        text(0, 1.1*max(ylim), ['Cell: ' num2str(iC)], 'fontweight', 'bold');
        pause(1);
    end
    
    [rate_map_B(:,:,iC)] = MS_decon_rate_map(ms.Binary(move_idx,iC), ms.time(move_idx), behav_aligned.position(move_idx,:), X_bins,Y_bins,0,smooth_sd);
    
end

%% same thing but split the data and check the corr

for iC = size(ms.RawTraces,2):-1:1
    [rate_map_1(:,:,iC)] = MS_decon_rate_map(deconv(move_idx_1,iC), ms.time(move_idx_1), behav_aligned.position(move_idx_1,:), X_bins,Y_bins, 0, smooth_sd);
    [rate_map_2(:,:,iC)] = MS_decon_rate_map(deconv(move_idx_2,iC), ms.time(move_idx_2), behav_aligned.position(move_idx_2,:), X_bins,Y_bins, 0, smooth_sd);

    [rate_map_B_1(:,:,iC)] = MS_decon_rate_map(ms.Binary(move_idx_1,iC), ms.time(move_idx_1), behav_aligned.position(move_idx_1,:), X_bins,Y_bins,0,smooth_sd);
    [rate_map_B_2(:,:,iC)] = MS_decon_rate_map(ms.Binary(move_idx_2,iC), ms.time(move_idx_2), behav_aligned.position(move_idx_2,:), X_bins,Y_bins,0,smooth_sd);

    this_map_1 = rate_map_1(:,:,iC);
    this_map_1(isnan(this_map_1)) = 0;
    this_map_2 = rate_map_2(:,:,iC);
    this_map_2(isnan(this_map_2)) = 0;
        split_corr(iC) = corr2(this_map_1, this_map_2); 

    this_map_1 = rate_map_B_1(:,:,iC);
    this_map_1(isnan(this_map_1)) = 0;
    this_map_2 = rate_map_B_2(:,:,iC);
    this_map_2(isnan(this_map_2)) = 0;
    
    split_corr_B(iC) = corr2(this_map_1, this_map_2); 

end


%% plot all cells rate map
% 
% x_tic = 0:3:size(occ_map,1)*3; x_tic = x_tic(1:end-1);
% y_tic = 0:3:size(occ_map,2)*3; y_tic = y_tic(1:end-1);
% %subplots
% N = 3; 
% M = 4;
% fig_range = 1:N*M:size(deconv,2)+1;
% 
% this_sub = 0;
% 
% for iC = 1:size(deconv,2) 
%     if ismember(iC, fig_range)
%         fig_n = find(iC == fig_range);
%         this_sub = 0;
%     end
%     this_sub = this_sub+1;
%     figure(350+fig_n)
%     
%     subplot(N,M,this_sub)
%     
%     imagesc(x_tic, y_tic,rate_map(:,:,iC)')
%     set(gca, 'YDir', 'normal');
%     set(gca, 'xtick', [], 'ytick', [])
% 
%         title({['Cell: ' num2str(iC)], ['corr: ' num2str(split_corr(iC),3)]},'color', [0.5 0.5 0.5])
% end

%% remove cells with less than .5 corr

split_thresh = 0.5; 
keep_idx = split_corr >split_thresh;
cell_ids = 1:length(split_corr);

fprintf('<strong>%d/%d</strong> cells passed the split half corr threshold (%0.2f)\n', length(keep_idx), length(split_corr), split_thresh)


% remove deconv with no spks
zero_keep_idx = zeros(1,size(rate_map,3));

for iC = size(rate_map,3):-1:1
    if sum(rate_map(:,:,iC),'all') ~=0
        zero_keep_idx(iC) = 1; 
    end
end
        cell_ids(~keep_idx | ~zero_keep_idx) = []; 

        rate_map(:,:,~keep_idx | ~zero_keep_idx) = [];
    
%% crate a shuffle set of maps for the null EMD distribution
Fs = mode(diff(ms.time));
nShuff  = 10000;
for iS = nShuff:-1:1
    
    dis_cells(iS) = floor(MS_randn_range(1,1,1,size(deconv,2)));  % track the cell selection
    dis_shuf(iS) = round(MS_randn_range(1,1, 4*Fs, length(deconv)+4*Fs)); % track the shift to ensure normal.
    
    this_shuff = circshift(deconv(:,dis_cells(iS)),dis_shuf(iS));  % shift the data, 4s minimum.
    
    [shuff_rate_map(:,:,iS)] = MS_decon_rate_map(this_shuff(move_idx), ms.time(move_idx), behav_aligned.position(move_idx,:), X_bins, Y_bins, 0, smooth_sd);
    
end

zero_keep_idx = zeros(1,size(shuff_rate_map,3));

for iC = size(shuff_rate_map,3):-1:1
    if sum(shuff_rate_map(:,:,iC),'all') ~=0
        zero_keep_idx(iC) = 1; 
    end
end

shuff_rate_map(:,:,~zero_keep_idx) = [];

% check of needed
% hist(dis_shuf,100)
% hist(dis_cells,100)

%% look at a specific cell


iC  =10;
if ishandle(102)
    close(102)
end
figure(iC)
set(gcf, 'position', [680 466 1085 505])
subplot(2,8,1:8)
hold on
plot(ms.time, ms.detrendRaw(:,iC), 'k');
plot(ms.time, deconv(:,iC) -.3,'color', [0.3467    0.5360    0.6907]);
plot(ms.time, denoise(:,iC), 'color', [ 0.9153    0.2816    0.2878 ] )
plot(ms.time, (ms.Binary(:,iC)*.1)-.41, 'color', [0.4416    0.7490    0.4322])
legend('Raw', 'Decon', 'Denoised', 'Binary', 'Orientation', 'horizontal','Location','North');
xlim([ms.time(1) ms.time(end)])

subplot(2,8, 9:10)
hold on
spk_idx = find(deconv(:,iC));
b_idx = find(ms.Binary(:,iC));

plot(behav_aligned.position(:,1), behav_aligned.position(:,2), 'color', [.8 .8 .8]);
plot(behav_aligned.position(b_idx,1), behav_aligned.position(b_idx,2),'.','markersize',20, 'color', [0.4416    0.7490    0.4322]);
plot(behav_aligned.position(spk_idx,1), behav_aligned.position(spk_idx,2),'.', 'color', [0.3467    0.5360    0.6907]);
xlim([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]);
ylim([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))]);
title(['Cell: ' num2str(iC)  'Slit corr: ' num2str(split_corr(iC),2)])

x_tic = 0:3:size(occ_map,1)*3; x_tic = x_tic(1:end-1);
y_tic = 0:3:size(occ_map,2)*3; y_tic = y_tic(1:end-1);

subplot(2,8,[11 12])
imagesc(x_tic, y_tic,occ_map(:,:,iC)')
set(gca, 'YDir', 'normal');
title('occupancy')

subplot(2,8,[13 14 ])
imagesc(x_tic, y_tic,rate_map(:,:,iC)')
set(gca, 'YDir', 'normal');
title('Decon rate')

subplot(2,8,[15 16 ])
imagesc(x_tic, y_tic,rate_map_B(:,:,iC)')
set(gca, 'YDir', 'normal');
title('Binary rate')


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
    templates.(temps{ii}) = imgaussfilt( templates.(temps{ii}),smooth_sd);
    
    % normalize
    templates.(temps{ii}) = templates.(temps{ii})./max(templates.(temps{ii}),[], 'all');
    
    subplot(1,length(temps),ii);
    imagesc(templates.(temps{ii}))
    xlabel(temps{ii})
end

% % square the subplots.  Ugly but at least it's sqaure
% axesHandles = findobj(get(401,'Children'), 'flat','Type','axes');
% % Set the axis property to square
% axis(axesHandles,'square')
%% try some methods on the elife data
% load('/home/ecarmichael/Downloads/data_share.mat')
addpath('/home/ecarmichael/Downloads/EMD_Yilmaz')


Temp_b = zeros(size(rate_map(:,:,1))); 
Temp_b(1:end,1) = 1; Temp_b(1:end,end) = 0;
Temp_b(1,1:end) = 0; Temp_b(end,1:end) = 0;
Temp_b = imgaussfilt(Temp_b,smooth_sd);

% normalize
Temp_b = Temp_b./max(Temp_b,[], 'all');

% EMD_score = datamat.rat167.scoreEMD(1); 

%% Yilmaz method to get null distribution
temp_b_gray = mat2gray(Temp_b); 
null_fval = [];
for iS = size(shuff_rate_map,3):-1:1
    
    [~, null_fval(iS)] = test(mat2gray(shuff_rate_map(:,:,iS)./max(shuff_rate_map(:,:,iS),[],'all')),temp_b_gray);
    
end
null_1prc = prctile(null_fval, 1);
%% try the Yilmaz method
fval = [];
for iC = size(rate_map,3):-1:1
%     map_1 = squeeze(datamat.rat167.spatialMap(ii,1,:,:)); 

    
[~, fval(iC)] = test(mat2gray(rate_map(:,:,iC)./ max(rate_map(:,:,iC),[],'all')),mat2gray(Temp_b));

% fprintf('cell: %d  |  fval = %0.2f  | emd = %0.2f',iC, fval, datamat.rat167.scoreEMD(1,ii))
end

% ranked version
[fval_s, s_idx] = sort(fval); 


% border cells by van Wijngaarden et al. 2020 method. 
b_cells = find(fval < null_1prc); 
non_b_cell =s_idx(~ismember(s_idx, b_cells));
non_b_fval = fval(non_b_cell);

%% plot the EMD distributions
% start with the null to get 1st percentile. 
figure(301)
subplot(3,5,6:10)
histogram(null_fval,25, 'facecolor', [0.8 0.8 0.8])
vline(null_1prc, '--k')
x_range= xlim; 
text(min(xlim), .8*max(ylim), ['1^s^t prc = ' num2str(null_1prc,3)])


% then do real data
[~, all_bins] = histcounts(fval, 25);
[count] = histcounts(fval(non_b_cell), all_bins);

subplot(3,5,1:5)
histogram('BinEdges', all_bins,'BinCounts', count, 'facecolor', [0.3467 0.5360 0.6907], 'facealpha', .2)
vline(null_1prc, '--k')
wide_xlim = [min([x_range, xlim]), max([x_range, xlim])];
xlim(wide_xlim)
hold on

subplot(3,5,1:5)
[count, ~] = histcounts(fval(b_cells),all_bins);
histogram('BinEdges', all_bins,'BinCounts', count, 'facecolor', [0.3467 0.5360 0.6907],'facealpha', 1)
xlim(wide_xlim)
text(min(xlim), .8*max(ylim), [num2str(length(b_cells)) '/' num2str(length(fval)) 'Border cells'])

 subplot(3,5,11)
    imagesc(x_tic, y_tic,Temp_b')
    set(gca, 'YDir', 'normal');
    set(gca, 'xtick', [], 'ytick', [])
    title('Template', 'color', [0.3467 0.5360 0.6907])

% add in some example cells
for iC = 1:length(b_cells)
    if iC > 2 % plot a max of 3 examples. 
        continue
    end
    subplot(3,5,11+iC)
    imagesc(x_tic, y_tic,rate_map(:,:,s_idx(iC))')
    set(gca, 'YDir', 'normal');
    set(gca, 'xtick', [], 'ytick', [])
    title({['Cell: ' num2str(s_idx(iC))], [' EMD = ' num2str(fval_s(iC),3)]}, 'color', [0.3467 0.5360 0.6907])

end

% plot a close non-border cell


 subplot(3,5,14)
 imagesc(x_tic, y_tic,rate_map(:,:,non_b_cell(1))')
 set(gca, 'YDir', 'normal');
 set(gca, 'xtick', [], 'ytick', [])
 title({['Cell: ' num2str(non_b_cell(1))], [' EMD = ' num2str(non_b_fval(1),3)]}, 'color', [.6 .6 .6])

 % and a middle of the pack non-border cell
 subplot(3,5,15)
 imagesc(x_tic, y_tic,rate_map(:,:,non_b_cell(floor(length(non_b_cell)/2)))')
 set(gca, 'YDir', 'normal');
 set(gca, 'xtick', [], 'ytick', [])
 title({['Cell: ' num2str(non_b_cell(floor(length(non_b_cell)/2)))], [' EMD = ' num2str(non_b_fval(floor(length(non_b_cell)/2)),3)]}, 'color', [.8 .8 .8])

if isfolder(save_dir)
 saveas(gcf, 'Border_summary.fig'); 
 saveas(gcf, 'Border_summary.png');
end
    

%% plot some example cells that pass the 1st prc of null

figure(333) % 
for iC = 1:length(b_cells)
    if iC > 3 % plot a max of 3 examples. 
        continue
    end
   subplot(3,5,iC)
   
    imagesc(x_tic, y_tic,rate_map(:,:,b_cells(iC))')
    set(gca, 'YDir', 'normal');
    title(['EMD = ' num2str(fval(b_cells(iC)),3)])
    
    
end

%% plot all cells with EMD score. 
%subplots
N = 3; 
M = 4;
fig_range = 1:N*M:size(deconv,2)+1;

this_sub = 0;



for iC = 1:size(deconv,2) 
    if ismember(iC, fig_range)
        fig_n = find(iC == fig_range);
        this_sub = 0;
    end
    this_sub = this_sub+1;
    figure(350+fig_n)
    
    subplot(N,M,this_sub)
    
    imagesc(x_tic, y_tic,rate_map(:,:,s_idx(iC))')
    set(gca, 'YDir', 'normal');
    set(gca, 'xtick', [], 'ytick', [])
    
    if fval_s(iC) < null_1prc
        set(gca, 'XColor', [0.3467 0.5360 0.6907],'YColor', [0.3467 0.5360 0.6907], 'linewidth', 3);
        title({['Cell: ' num2str(s_idx(iC))], [' EMD = ' num2str(fval_s(iC),3)]}, 'color', [0.3467 0.5360 0.6907])
    else
        title({['Cell: ' num2str(s_idx(iC))] ,[' EMD = ' num2str(fval_s(iC),3)]},'color', [0.5 0.5 0.5])
    end
end


%% try the EMD
im1 = map_1; 
im2 = Temp_b;

R= size(im1,1);
C= size(im1,2);
if (~(size(im2,1)==R&&size(im2,2)==C))
    error('Size of images should be the same');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMD input
% Each unit of gray-level (between 0 and 255) is a unit of mass, and the
% ground distance is a thresholded distance. This is similar to:
%  A Unified Approach to the Change of Resolution: Space and Gray-Level
%  S. Peleg and M. Werman and H. Rom
%  PAMI 11, 739-742
% The difference is that the images do not have the same total mass 
% and we chose a thresholded ground distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COST_MULT_FACTOR= 1000;
THRESHOLD= 3*COST_MULT_FACTOR;
D= zeros(R*C,R*C,'int32');
j= 0;
for c1=1:C
    for r1=1:R
        j= j+1;
        i= 0;
        for c2=1:C
            for r2=1:R
                i= i+1;
                D(i,j)= min( [THRESHOLD (COST_MULT_FACTOR*sqrt((r1-r2)^2+(c1-c2)^2))] );
            end
        end
    end
end
extra_mass_penalty= int32(0);
flowType= int32(3);

P= int32(im1(:));
Q= int32(im2(:));

[emd_hat_gd_metric_mex_val, score] = emd_hat_gd_metric_mex(P,Q,D,extra_mass_penalty);


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

