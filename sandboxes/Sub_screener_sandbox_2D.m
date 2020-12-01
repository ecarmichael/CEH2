
%% Sub screening 2d sandbox

close all
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;

if ismac
    %     PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    %     PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    %     PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    %     PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    %     PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    %
elseif strcmp(os, 'GLNXA64')
    
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.data_dir = '/mnt/Data/Williams_Lab/II_classification/msbehavplace/ck2cre1'; % where to find the raw data
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    %         PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/'; % raw data location.
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/II_classification/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/II_classification/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    %
else
    PARAMS.data_dir = 'J:\Williams_Lab\II_classification'; % where to find the raw data
    PARAMS.raw_data_dir = 'J:\Williams_Lab\II_classification'; % raw data location.
    PARAMS.inter_dir = 'J:\Williams_Lab\II_classification\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\II_classification\Inter\Stats'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.code_seqnmf_dir = 'C:\Users\ecarm\Documents\GitHub\seqNMF'; % where the multisite repo can be found
    
end

% colours
PARAMS.L_grey = [0.8 0.8 0.8];
PARAMS.D_grey = [0.2 0.2 0.2];
PARAMS.blue = [0.3639    0.5755    0.7484];
PARAMS.red = [0.9153    0.2816    0.2878];
PARAMS.green= [0.4416    0.7490    0.4322];
PARAMS.gold = [1.0000    0.5984    0.2000];

rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

% make the Inter and Stats dir if they don't exist.
if ~exist(PARAMS.stats_dir,'dir')
    mkdir(PARAMS.stats_dir)
end

% % try the newer NLX loaders for UNIX
% [~, d] = version;
% if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
%     rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
%     addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
%     disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
%     which Nlx2MatCSC
% end

clear d os

%%  load data

load('behav.mat')
load('ms.mat')
ms = MS_msExtractBinary_detrendTraces(ms, 2);

% cfg_rm.remove_idx = [17 26];
% ms = MS_Remove_trace(cfg_rm, ms);



figure(101)
cfg_plot = [];
cfg_plot.view =[0 75];
cfg_plot.plot_type = '2d';
% cfg_plot.colors = parula(size(ms.Binary,2));
MS_plot_ca(cfg_plot, ms)
xlabel('time (s)')



%% align behaviour and Ca

if behav.time(end) ~= ms.time(end) || length(behav.time) ~= length(ms.time)
    fprintf('<strong> %s </strong>: behaviour and Ca are not the same length ro end time.  attempting alignment \n', mfilename);
    
    
    behav_aligned = MS_align_data(behav, ms);
    %     behav_aligned = MS_align_data(behav, ms);
    
end

%smooth speed
behav_aligned.speed = smooth(behav_aligned.speed, 3*mode(diff(ms.time)));

cfg.min_speed = 5; % in cm/s
movement_idx = behav_aligned.speed >cfg.min_speed; % get times when the animal was moving.





%% plot basics for each cell
for iC = 1:size(ms.Binary,2) % loop through cells
pthreshold = 0.05; % value for pvalue cut off;
M = 3; % rows
N = 5; % columns
%   fig{iC} = figure('Visible', 'off');

%get binary 'event times' to be plotted as dots
t_binary = find(ms.Binary(:,iC) ==1);

%%% title information
subplot(M, N, [N-2 (N*2)-2]) % title information. Top right corner,
text(0, 10, ['Cell id: ' num2str(iC)])
text(0, 8, ['Binary thresh: ' num2str(ms.Binary_threshold)])
ylim([0 10])
axis off


%%% raw trace
subplot(M, N, 1:3)
plot(ms.time/1000, ms.RawTraces(:,iC), 'color', PARAMS.blue)
xlim([ms.time(1)/1000 ms.time(end)/1000]);
xlabel('time(s)');
ylabel('dF/F');
hline(mean(ms.RawTraces(:,iC))+2*std(ms.RawTraces(:,iC)));



%%% X Y position
subplot(M, N, N+1:N+3)
hold on
plot(behav_aligned.time/1000, behav_aligned.position(:,1), 'color', PARAMS.L_grey)
plot(behav_aligned.time/1000, behav_aligned.position(:,2), 'color', PARAMS.D_grey)

%%% update position in time with binary 'spikes'
plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red);
plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,2),'.', 'color', PARAMS.red);

% plot(behav_aligned.time/1000, behav_aligned.position(:,2),'color', PARAMS.blue)
ylabel('linear position')
xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
% legend({'x', 'y'})


%%% speed info
subplot(M, N, N*2+1:N*2+3)
hold on
plot(behav_aligned.time/1000, behav_aligned.speed, 'color', PARAMS.L_grey)
plot(behav_aligned.time(movement_idx)/1000, behav_aligned.speed(movement_idx),'.', 'color', PARAMS.gold, 'markersize', 1)
% legend('Speed', 'box', 'off')

xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
ylabel('Speed cm/s')
xlabel('time (s)')

%%% update speed in time with binary 'spikes'
subplot(M, N, N*2+1:N*2+3)
plot(behav_aligned.time(t_binary)/1000,behav_aligned.speed(t_binary,1),'.', 'color', PARAMS.red);



% %%% orientation info
    % subplot(M, N, N*3+1:N*3+3)
% plot(behav_aligned.time/1000,ones(size(behav_aligned.time)), 'color', 'w')
% hold on
% text(behav_aligned.time(floor(length(behav_aligned.time)/3))/1000, pi, 'HD placeholder')
%     % ylabel('HD')
% ylim([-pi pi])
% set(gca, 'ytick', [-pi pi], 'yticklabel', {'-pi' 'pi'})

%%% plot the binary times on the position
subplot(M, N, [N+4]) % N*4+4:N*4+6
hold on
plot(behav_aligned.position(:,2), behav_aligned.position(:,1), 'color', PARAMS.L_grey)
xlim([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))])
ylim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]));

% put dots on positions when the cell was active.
plot(behav_aligned.position(t_binary,2), behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red)
xlabel('position (cm)');
ylabel('position (cm)');
ylim(round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]))
set(gca, 'ytick', round([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))]));
% set(gca, 'yticklabel', num2str(roundn([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))],2)))
% get the transient/position values
% tran_x = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,1),ms.time(t_binary),'linear');
% tran_y = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,2),ms.time(t_binary),'linear');
%
% plot(tran_x,tran_y,'.', 'color', PARAMS.red);




%%% add the SPF for this cell.
subplot(M, N, N) % spf with centroid.
Spr = winter(32);
colormap([0 0 0 ; Spr(16:end,:)]);
% c_lim = [0.2*max(max(ms.PeakToNoiseProj)), max(max(ms.PeakToNoiseProj))]; % helps clean up the projection by increasing the floor of the colormap to take in the lowest x% of the data
% imagesc(ms.PeakToNoiseProj, c_lim)
MS_plot_all_SFPs(flipdim(ms.SFPs,3)); % custom function to plot all the SFPs on top of each other.  Cleaner than ms.PeakToNoiseProj.
hold on
quiver(ms.Centroids(iC,2)-30,ms.Centroids(iC,1)+5, 22,-3,-10,'color', 'w', 'linewidth', 2, 'MaxHeadSize', 5); % add an arrow pointing to the current cell.
%     scatter(ms.Centroids(iC,2), ms.Centroids(iC,1),60,'w', 'o','LineWidth',.5); % put a circle around the current cell.


%% add in the 2D place/spatial information?
    subplot(M, N, N+5) % N*4+4:N*4+6
    cfg.bin_size = 3;
    X_bins = min(behav_aligned.position(:,1)):cfg.bin_size:max(behav_aligned.position(:,1));
    X_bin_centers = X_bins +  cfg.bin_size/2;
    X_bin_centers = X_bin_centers(1:end-1);
    % same for Y bins
    Y_bins = min(behav_aligned.position(:,1)):cfg.bin_size:max(behav_aligned.position(:,2));
    Y_bin_centers = Y_bins +  cfg.bin_size/2;
    Y_bin_centers = Y_bin_centers(1:end-1);
    
    [MI, posterior, occupancy, p_active, likelihood] = MS_get_spatial_information_2D(ms.Binary(:,iC),behav_aligned.position, X_bins, Y_bins );
    % GE method for confirmation
%     [MI, PDF, occupancy_map, prob_being_active, tuning_map] = extract_2D_information(ms.Binary(:,iC),behav_aligned.position, X_bins, Y_bins, 1:length(behav_aligned.position));

    imagesc(X_bin_centers,Y_bin_centers,likelihood);
    axis xy    
    xlabel('position (cm)');
    ylabel('position (cm)');

    
    %%
% Bootstramp method
nShuff = 1000;
actual_bootstrap_tuning_curve = zeros(length(X_bin_centers),length(Y_bin_centers), nShuff);
shuffled_bootstrap_tuning_curve = zeros(length(bin_centers), nShuff);
for iShuff = nShuff:-1:1
    split_ts = ceil(MS_randn_range(1,1,1,length(ms.time)));
    this_shuff = [ms.Binary(end-split_ts+1:end,iC); ms.Binary(1:end-split_ts,iC)]; % cut the data at a point and put the ends together.
    
    bootstrap_ts = right_idx;
    for ii = 1:length(bootstrap_ts)
        if bootstrap_ts(ii) == 1 && rand < 0.5
            bootstrap_ts(ii) = 0;
        end
    end
    
    % Compute the actual tuning curve using a bootstrapped sample
    %         [actual_bootstrap_MI(iShuff), actual_bootstrap_PDF(:,iShuff), ~, actual_bootstrap_prob_being_active(iShuff), actual_bootstrap_tuning_curve(:,iShuff) ] = MS_get_spatial_information(ms.Binary(:,iC), behav_aligned.position(:,1), bins, split_ts);
    [actual_bootstrap_MI(iShuff), actual_bootstrap_PDF(:,iShuff), ~, actual_bootstrap_prob_being_active(iShuff), actual_bootstrap_tuning_curve(:,iShuff)] = MS_get_spatial_information(ms.Binary(bootstrap_ts,iC),ms.time(bootstrap_ts), behav_aligned.position(bootstrap_ts,1), bins);
    
    % Compute the shuffled tuning curve using the same bootstrapped sample
    [shuffled_bootstrap_MI(iShuff), shuffled_bootstrap_PDF(:,iShuff), ~, shuffled_bootstrap_prob_being_active(iShuff), shuffled_bootstrap_tuning_curve(:,iShuff)] = MS_get_spatial_information(this_shuff(bootstrap_ts), ms.time(bootstrap_ts),behav_aligned.position(bootstrap_ts,1), bins);
end
%     hold on
%     plot(bin_centers,shuffled_bootstrap_tuning_curve, 'r')
%     plot(bin_centers,actual_bootstrap_tuning_curve, 'k')


% Find the 95% confidence interval
sorted_BS_tuning_curves = sort(actual_bootstrap_tuning_curve,2);
CI_idx_loc = 0.95*nShuff/2;
median_idx = round(nShuff/2);
upper_CI95_idx = median_idx+CI_idx_loc;
lower_CI95_idx = median_idx-CI_idx_loc;

% This will make sure that upper and lower bounds are withing the actual bootstrap data
upper_CI95_idx(upper_CI95_idx > nShuff) = nShuff;
upper_CI95_idx(upper_CI95_idx < 1) = 1;
lower_CI95_idx(lower_CI95_idx > nShuff) = nShuff;
lower_CI95_idx(lower_CI95_idx < 1) = 1;

upper_CI95 = sorted_BS_tuning_curves(:,upper_CI95_idx);
lower_CI95 = sorted_BS_tuning_curves(:,lower_CI95_idx);

plot(bin_centers, actual_bootstrap_tuning_curve, 'color', [0.8 0.8 0.8], 'Linewidth', 0.5)

plot(bin_centers,likelihood, 'k', 'Linewidth', 1)
plot(bin_centers,upper_CI95, 'r', 'Linewidth', 1)
plot(bin_centers,lower_CI95, 'r', 'Linewidth', 1)
ylabel('p(A|S) R')

%%
%     % add in shuffle (based on Guillaume's CaImDecoding)
%     subplot(M, N, N*4+1:N*4+3)
%     nShuff = 1000;
%     for iSuff = nShuff:-1:1
%         split_ts = ceil(MS_randn_range(1,1,1,length(ms.time)));
%         this_shuff = [ms.Binary(end-split_ts+1:end,iC); ms.Binary(1:end-split_ts,iC)]; % cut the data at a point and put the ends together.
%         [~, ~, ~, ~, shuff_TC(:, iSuff)] = MS_get_spatial_information(this_shuff, ms.time, behav_aligned.position(:,1), bins);
%
%     end
%     pval = sum(shuff_TC> likelihood,2)/nShuff; %  p-value, supra-threshold tests
%     sig_vals = likelihood;
%     sig_vals(pval > pthreshold) = 0;
%     hold on
%     plot(bin_centers, likelihood, 'k');
%     plot(bin_centers, sig_vals, 'r');
%     ylabel('p val')




%% add speed mod score
    subplot(M, N, N*2+4) % N*4+4:N*4+6

% make speed bins
keep_idx = behav_aligned.speed >= cfg.min_speed; 
cfg.speed_bin_size = 2;
 S_bins =  cfg.min_speed:cfg.speed_bin_size:max(behav_aligned.speed(keep_idx,1));
    S_bin_centers = S_bins +  cfg.bin_size/2;
    S_bin_centers = S_bin_centers(1:end-1);
    
    speed_hist = hist(behav_aligned.speed(keep_idx), S_bin_centers);
    speed_active_hist = hist(behav_aligned.speed(t_binary), S_bin_centers); 
    yyaxis left
    bar(S_bin_centers, speed_hist./mode(diff(behav_aligned.time)), 'facecolor', PARAMS.D_grey, 'edgecolor', PARAMS.D_grey)
        ylabel('time at speed (s)')

%         plot(S_bin_centers, speed_active_hist./mode(diff(behav_aligned.time)), 'color', PARAMS.red)
        yyaxis right
    plot(S_bin_centers, speed_active_hist./speed_hist, 'color', PARAMS.red)
%     legend({'occupancy', 'acitve'});
    xlabel('speed (cm/s)')
    ylabel('p active')
    
            [S_MI, ~, ~,~, S_likelihood] = MS_get_spatial_information(ms.Binary(:,iC),ms.time, behav_aligned.speed, S_bins);

    
    % get boostrapped values for speed using MI
    nShuff = 1000;
    actual_bootstrap_tuning_curve = zeros(length(S_bin_centers), nShuff);
    shuffled_bootstrap_tuning_curve = zeros(length(S_bin_centers), nShuff);
    for iShuff = nShuff:-1:1
        split_ts = ceil(MS_randn_range(1,1,1,length(ms.time)));
        this_shuff = [ms.Binary(end-split_ts+1:end,iC); ms.Binary(1:end-split_ts,iC)]; % cut the data at a point and put the ends together.
        
        % get a set of indicies to include for this shuffle
        bootstrap_ts = ones(1,length(this_shuff));
        for ii = 1:length(bootstrap_ts)
            if bootstrap_ts(ii) == 1 && rand < 0.5
                bootstrap_ts(ii) = 0;
            end
        end
        bootstrap_ts = logical(bootstrap_ts); 
        % Compute the actual tuning curve using a bootstrapped sample
        [actual_bootstrap_MI(iShuff), actual_bootstrap_PDF(:,iShuff), ~, actual_bootstrap_prob_being_active(iShuff), actual_bootstrap_tuning_curve(:,iShuff)] = MS_get_spatial_information(ms.Binary(bootstrap_ts,iC),ms.time(bootstrap_ts), behav_aligned.speed(bootstrap_ts), S_bins);
        
        % Compute the shuffled tuning curve using the same bootstrapped sample
        [shuffled_bootstrap_MI(iShuff), shuffled_bootstrap_PDF(:,iShuff), ~, shuffled_bootstrap_prob_being_active(iShuff), shuffled_bootstrap_tuning_curve(:,iShuff)] = MS_get_spatial_information(this_shuff(bootstrap_ts), ms.time(bootstrap_ts),behav_aligned.speed(bootstrap_ts), S_bins);
    end
%     hold on
%     plot(bin_centers,shuffled_bootstrap_tuning_curve, 'r')
%     plot(bin_centers,actual_bootstrap_tuning_curve, 'k')


% Find the 95% confidence interval
sorted_BS_tuning_curves = sort(actual_bootstrap_tuning_curve,2);
CI_idx_loc = 0.95*nShuff/2;
median_idx = round(nShuff/2);
upper_CI95_idx = median_idx+CI_idx_loc;
lower_CI95_idx = median_idx-CI_idx_loc;

% This will make sure that upper and lower bounds are withing the actual bootstrap data
upper_CI95_idx(upper_CI95_idx >  nShuff) =  nShuff;
upper_CI95_idx(upper_CI95_idx < 1) = 1;
lower_CI95_idx(lower_CI95_idx >  nShuff) =  nShuff;
lower_CI95_idx(lower_CI95_idx < 1) = 1;

upper_CI95 = sorted_BS_tuning_curves(:,upper_CI95_idx);
lower_CI95 = sorted_BS_tuning_curves(:,lower_CI95_idx);

plot(S_bin_centers, actual_bootstrap_tuning_curve', 'color', [0.8 0.8 0.8], 'Linewidth', 0.5)
hold on
plot(S_bin_centers,S_likelihood, 'color', PARAMS.blue, 'Linewidth', 2)
plot(S_bin_centers,upper_CI95, 'color', PARAMS.red, 'Linewidth', 2)
plot(S_bin_centers,lower_CI95, 'color', PARAMS.red, 'Linewidth', 2)

%     [AX,H1,H2] =plotyy(S_bin_centers, speed_hist,S_bin_centers, speed_active_hist, 'bar', 'bar');
% set(H1,'FaceColor','r') % a
% set(H2,'FaceColor','b') % b

    % use 1d extract based on speed bins. 
    
    
%% customize figure stuff

pos = get(gcf, 'position');
set(gcf, 'position', [pos(1)-pos(1)*.6 pos(2)-pos(2)*.6 pos(3)*1.8 pos(4) *1.6])


% pause(3)
%     close(100)
[full,this_dir]=fileparts(pwd);
[~,this_parent] = fileparts(full);
mkdir([PARAMS.inter_dir  this_parent filesep this_dir]);
saveas(gcf, [PARAMS.inter_dir  this_parent filesep this_dir filesep 'Cell_' num2str(iC) '_Spatial_info.fig'])
saveas(gcf, [PARAMS.inter_dir  this_parent filesep this_dir filesep 'Cell_' num2str(iC) '_Spatial_info.png'])


close(iC)
end


%% make a plot of population level activity 

%population speed modulation

% bin pop binary signal into .5s bins (Chen 2015) to get % of cells active
% at a given time. get Peasron cor coeff between speed and pop activity.
% they normalized by mean and SD for both the speed vec and pop activity
% vec. 








