function [rem_z, Ass_z, time_proj_pos_place, wake_time_proj_rem] = sandbox_PCA_ICA_no_fig(fname)
% sandbox_PCA/ICA


if strcmp(computer, 'GLNXA64')
    
    codebase_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/home/williamslab/Documents/Github/CEH2';
    oasis_dir = '/home/williamslab/Documents/Github/OASIS_matlab';
    
    code_dir = '/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '/home/williamslab/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = '/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter';
    rem_dir = '/home/williamslab/Williams Lab Dropbox/Eric Carmichael/JisooProject2020/2020_Results_aftercutting/Across_episodes/Inter';
    decode_dir = [data_dir filesep 'decoding'];
    
    this_process_dir ='/home/williamslab/Williams Lab Dropbox/Eric Carmichael/JisooProject2020/2020_Results_aftercutting/';
    
    
else
    
    codebase_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
    ca_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
    
    code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';
    
    RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
    rem_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter';
    decode_dir = [data_dir filesep 'decoding'];
    this_process_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\JisooProject2020\2020_Results_aftercutting';
    
end

restoredefaultpath

cd(oasis_dir)
addpath(genpath(oasis_dir));
oasis_setup

addpath(genpath(ca_dir));
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(code_dir)

cd(data_dir)

rng(123, 'twister')

%%

% load('ms_trk.mat')
% load('behav_DLC.mat')
this_sess = fname;

load(this_sess)

behav = MS_align_data(behav,ms);

move_idx = behav.speed > 2.5;

%% check for place cell metrics
dir_parts = strsplit(this_sess, filesep);
parts = strsplit(dir_parts{end}, '_');
task = parts{2};
if contains(task, 'HATD6')
    task = 'HATDSwitch';
end
subject = parts{1};


if exist([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'], 'file')
    load([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'])
elseif exist([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'], 'file')
    load([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'])
else
    rem_z = NaN; Ass_z = NaN;  time_proj_pos_place = NaN; wake_time_proj_rem = NaN;
    return
end
%%

ms_trk = ms;
keep_idx = zeros(1,size(ms_trk.RawTraces,2));
keep_idx(1:floor(size(ms_trk.RawTraces,2)*.66)) = 1;

remove_cell_id = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'RawTraces';
ms_trk_rem = MS_Remove_trace(cfg_rem, ms_trk);


ms_trk_rem = MS_append_deconv(ms_trk_rem, 1);

% remove inactive cells

keep_idx = sum(ms_trk_rem.deconv, 1) >0;

remove_cell_id_decon = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv';
ms_trk_rem = MS_Remove_trace(cfg_rem, ms_trk_rem);

%% grab the spatial tuning properties.

% if contains(task, 'HATD')% workaround for naming of HAT + D
%     load([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'])
% %     load([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'SA.mat'])
% else
%     load([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'])
% end

place = [];
for iC = length(spatial_analysis.bin):-1:1
    if iscell(spatial_analysis.bin{iC,1}.PlaceFieldCentroid)
        temp = cell2mat(spatial_analysis.bin{iC,1}.PlaceFieldCentroid);
    else
        temp = spatial_analysis.bin{iC,1}.PlaceFieldCentroid;
    end
    place.centroids(iC) = temp(1);
    
    % is it a place cell?
    place.is(iC) = spatial_analysis.bin{iC,1}.IsPlaceCell;
    
    place.map(iC,:) = mean(spatial_analysis.bin{iC,1}.PlaceField,1)/max(mean(spatial_analysis.bin{iC,1}.PlaceField,1)); % get the 1d place field and normalize.
end

place.centroids(remove_cell_id)= [];
place.centroids(remove_cell_id_decon)= [];

place.is(remove_cell_id)= [];
place.is(remove_cell_id_decon)= [];

place.map(remove_cell_id,:)= [];
place.map(remove_cell_id_decon,:)= [];

bin = 80/size(place.map,2);
p_bins = 10:bin:90;
p_bins = p_bins(1:end-1)+bin/2;
% see if there are any anxiety cells


%% follow grosmark et al. method of deconv preprocessing
Csp = ms_trk_rem.deconv./ms_trk_rem.denoise;
Csp = Csp > 0.01;
ms_trk_rem.Csp = Csp;

% cfg_plot.Ca_type = 'RawTraces';
% cfg_plot.plot_type = '2d';
% MS_plot_ca(cfg_plot, ms_trk_rem)

%% bin and convolve
binsize = 0.1; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window

gau_z = zscore(gau_sdf, [], 2);
%% plot the gaussian smoothed Csp
% figure(101)
% ax(1)= subplot(5,1,1);
% hold on
% plot(ms_trk.time, behav.position(:,1));
% plot(ms_trk.time, behav.position(:,2));
% xlim([ms_trk.time(1) ms_trk.time(end)])
%
% ax(2) = subplot(5,1,2:5);
% imagesc(ms_trk.time, 1:size(Csp,2), gau_z')
%
% linkaxes(ax, 'x')


% %% optional bin the time into something larger than a frame.


% bin_s = 1; % time in seconds
%
% tvec = round(ms.time(1))/1000:bin_s:round(ms.time(end)/1000);
%
%
% data_int = []; data_sdf = [];
% for ii = size(gau_z,2):-1:1
%     data_int(:,ii) = interp1(ms.time/1000, gau_z(:,ii), tvec);
%     data_sdf (:,ii) = interp1(ms.time/1000, gau_sdf(:,ii), tvec);
% end
%
%
% % spk_count = histcounts(



%% try again with histc
binsize = .5;
tbin_edges = ms.time(1)/1000:binsize:ms.time(end)/1000; % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_h = [];
for ii = size(Csp,2):-1:1
    
    this_cell = ms.time(find(Csp(: ,ii) & move_idx))/1000;
    
    spk_count = histc(this_cell,tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    
    data_h(:,ii) = spk_count;
end



tvec = tbin_centers;


%% convert data if needed.
% gau_z_t = gau_z;
% gau_sdf_t = gau_sdf;
%
% gau_z = data_int;
% gau_sdf = data_int;

%% try the assembly code....
rng(123, 'twister')
Ass_Temp = assembly_patterns(data_h');

rng(123, 'twister')
time_proj = assembly_activity(Ass_Temp,data_h');


%% shuffle distribution for assemblies
rng(123,'twister')
nShuff = 200;

Ass_shuff = NaN(1,nShuff);
parfor iS = 1:nShuff
    tic
    shuff_data = NaN(size(data_h));
    for ic = 1:size(data_h,2)
        shuff_data(:,ic) = circshift(data_h(:,ic), floor(MS_randn_range(1,1,1,size(data_h,1))));
    end
    
    this_ass = assembly_patterns(shuff_data');
    
    
    %     for ii = size(this_ass,2):-1:1
    
    if sum(max(this_ass) > 0.2) >0
        Ass_shuff(iS) = sum(max(this_ass) > 0.2);
    else
        Ass_shuff(iS) = 0;
    end
    %     end
    fprintf('Shuff # %.0f found %.0f assemblies and took %2.2f seconds\n', iS, size(this_ass,2), toc)
end


% %% color code the assemblies
% Ass_sort = [];
% %figure(303);clf;  hold on
% for ii = size(Ass_Temp,2):-1:1
%
%     [Ass_sort(:,ii), this_idx] = sort(Ass_Temp(:,ii), 'descend');
%     thresh = prctile(Ass_Temp(:,ii), 95);
%     c_idx = Ass_Temp(this_idx,ii) > thresh;
%
%     %     plot(Ass_Temp(c_idx, ii))
%     %             stem(Ass_Temp(c_idx, ii))
%
% end

%% keep assemblies with strong activations
keep_idx = zeros(1,size(Ass_Temp,2));
Ass_pos_cells = cell(size(keep_idx));
for ii = size(Ass_Temp,2):-1:1
    
    if max(Ass_Temp(:, ii)) > 0.2
        keep_idx(ii) = 1;
        
        z_weight = zscore(Ass_Temp(:,ii));
        
        Ass_pos_cells{ii} = find(z_weight > 2);
        
    end
    
    
end

Ass_pos = Ass_Temp;
Ass_pos(:,~keep_idx) = [];
Ass_pos_cells(~keep_idx) = [];
Ass_idx = find(keep_idx);

time_proj_pos = time_proj;
time_proj_pos(~keep_idx,:) = [];

Ass_z = (size(Ass_pos,2) - mean(Ass_shuff))/std(Ass_shuff);
fprintf('%.0f Positive Assemblies detected. Chance level is %.1f. zscore = %.1fSD\n', size(Ass_pos,2), mean(Ass_shuff), Ass_z)
%% stem plot ensembles
% figure(301);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% figure(302);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

Ass_map = cell(size(Ass_pos,2),1);
Ass_mean_map = [];
Ass_pcells = Ass_map;
c_ord = parula(size(Ass_pos,2)+2);
%
% cmap = parula(256);
% cmap(1,:) = 0;
for ii = size(Ass_pos,2):-1:1
    %     figure(301)
    %     subplot(4, ceil(size(Ass_pos,2)/4),ii)
    %     hold on
    %     stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
    %     view(90,90)
    
    %     stem(Ass_pos_cells{ii}, Ass_pos(Ass_pos_cells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
    %         title(['Assembly #' num2str(ii)])
    
    
    %     figure(302)
    %         subplot(4, ceil(size(Ass_pos,2)/4),ii)
    %     hold on
    this_ass_map = []; these_place = zeros(size(Ass_pos_cells{ii}));
    fprintf('Assembly # %.0f had %.0f place cells\n', ii, sum(place.is(Ass_pos_cells{ii})))
    for jj = 1:length(Ass_pos_cells{ii})
        
        if place.is(Ass_pos_cells{ii}(jj))
            these_place(jj) = 1;
            place_int = interp1(p_bins,place.map(Ass_pos_cells{ii}(jj),:),  p_bins(1):1:p_bins(end));
            
            this_ass_map = [this_ass_map ; place_int];
        end
    end
    
    if ~isempty(this_ass_map)
        %         imagesc(p_bins(1):1:p_bins(end), 1:length(sum(these_place)),  this_ass_map)
        %         set(gca,'ytick',1:size(this_ass_map,1), 'YTickLabel',  Ass_pos_cells{ii}(logical(these_place))');
        %         xlim([p_bins(1) p_bins(end)])
        %         ylim([.5 size(this_ass_map,1)+.5])
        %         colormap(cmap)
        Ass_map{ii} = this_ass_map;
        Ass_pcells{ii} = Ass_pos_cells{ii}(logical(these_place));
        Ass_mean_map(ii,:) = mean(this_ass_map,1);
    else
        
        Ass_map{ii} = NaN(size(p_bins(1):1:p_bins(end)));
        Ass_pcells{ii} = NaN;
        Ass_mean_map(ii,:) = NaN(size(p_bins(1):1:p_bins(end)));
    end
    %     title(['Assembly #' num2str(ii)])
    
    
    
    %     ms_t = MS_append_sharp_SFPs(ms_trk);
    
    
    %     MS_plot_all_SFPs(imgaussfilt(ms_t.SFPs_sharp(:,:, Ass_pos_cells{ii}),2))
    
end
nan_idx = isnan(max(Ass_mean_map, [], 2))
Ass_idx(nan_idx) = [];
Ass_pos(:,nan_idx) = [];
Ass_mean_map(nan_idx,:) = [];
Ass_map(nan_idx) = [];
Ass_pcells(nan_idx) = [];

Ass_pos_cells_place = Ass_pos_cells;
time_proj_pos_place = time_proj_pos;

Ass_pos_cells_place(nan_idx) = [];
time_proj_pos_place(nan_idx,:) = [];


% Ass_pos = Ass_Temp;
% Ass_pos(:,~keep_idx) = [];
% Ass_pos_cells(~keep_idx) = [];
%
% time_proj_pos = time_proj;
% time_proj_pos(~keep_idx,:) = [];

%% sort the maps based on peak spatial representation
[~, idx] = max(Ass_mean_map, [], 2);
p_bins_int = p_bins(1):1:p_bins(end);

Ass_map_peak = p_bins_int(idx);
[Ass_map_peak, s_idx] = sort(Ass_map_peak);

Ass_idx = Ass_idx(s_idx);
Ass_pos = Ass_pos(:,s_idx);
Ass_mean_map = Ass_mean_map(s_idx,:);
Ass_map = Ass_map(s_idx);
Ass_pcells = Ass_pcells(s_idx);
Ass_pos_cells_place = Ass_pos_cells_place(s_idx);
time_proj_pos_place = time_proj_pos_place(s_idx,:);



rm_idx =  find(cellfun(@length, Ass_pcells) < 4);
Ass_map_peak(rm_idx) = [];
Ass_idx(rm_idx) = [];
Ass_pos(:,rm_idx) = [];
Ass_mean_map(rm_idx,:) = [];
Ass_map(rm_idx) = [];
Ass_pcells(rm_idx) = [];
Ass_pos_cells_place(rm_idx) = [];
time_proj_pos_place(rm_idx,:) = [];



%% get the number of significant reactivations during wake for


c_ord = cool(size(time_proj_pos_place,1)+4);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grab the REM data
sess = strsplit(this_sess, '_');
sub = sess{1};
sess = sess{2};
cd(rem_dir )

if strcmpi(sess, 'HATDSwitch')
    sess = 'HATSwitch';
end

sub_list = dir('PV*');
keep_idx = [];
for ii = 1:length(sub_list)
    if strfind(lower(sub_list(ii).name), lower(sub))
        keep_idx(ii) = 1;
    else
        keep_idx(ii) = 0;
    end
end

cd([rem_dir filesep sub_list(find(keep_idx)).name])

sess_list = dir('*T*');
keep_idx = [];
for ii = 1:length(sess_list)
    if strfind(lower(sess_list(ii).name), lower(sess))
        keep_idx(ii) = 1;
    else
        keep_idx(ii) = 0;
    end
end

cd(sess_list(find(keep_idx)).name)

if exist('all_detrendRaw_post_REM.mat', 'file')
    load('all_detrendRaw_post_REM.mat')
    load('all_RawTraces_post_REM.mat')
    load('all_binary_post_REM.mat')
else
    rem_z = NaN; wake_time_proj_rem = NaN;
    return
end


%remove cells that were excluded in the awake set.

all_detrendRaw_post_REM(:, remove_cell_id) = [];
all_detrendRaw_post_REM(:, remove_cell_id_decon) = [];

all_RawTraces_post_REM(:, remove_cell_id) = [];
all_RawTraces_post_REM(:, remove_cell_id_decon) = [];

all_binary_post_REM(:, remove_cell_id) = [];
all_binary_post_REM(:, remove_cell_id_decon) = [];
% load the decoding as well





%% deconvolve
addpath(genpath(oasis_dir))
fprintf('\n<strong>%s</strong>: deconvolving traces...\n', mfilename)
for iChan = size(all_detrendRaw_post_REM,2):-1:1
    tic;
    [denoise,deconv] = deconvolveCa(all_detrendRaw_post_REM(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
    toc;
    all_denoise(:,iChan) = denoise;    all_deconv(:,iChan) = deconv;
end
rmpath(genpath(oasis_dir))


% follow grosmark et al. method of deconv preprocessing
Csp_rem = all_deconv./all_denoise;
Csp_rem = Csp_rem > 0.01;


%%
rem_time =  0:1/mode(diff(ms.time)):(length(all_deconv)/mode(diff(ms.time)));
rem_time = rem_time(1:end-1);

binsize = .5;
tbin_edges_rem = 0:binsize:(length(all_deconv)/mode(diff(ms.time))); % vector of time bin edges (for histogram)
tbin_centers_rem = tbin_edges_rem(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_h_rem = [];
for ii = size(Csp_rem,2):-1:1
    
    this_cell = rem_time(find(Csp_rem(:,ii)));
    
    spk_count = histc(this_cell,tbin_edges_rem); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    
    data_h_rem(:,ii) = spk_count;
end



tvec_rem = tbin_centers_rem;

%% try the assembly code
rng(123, 'twister')
Ass_Temp_rem = assembly_patterns(data_h_rem');

rng(123, 'twister')
%using wake assemblies.
wake_time_proj_rem = assembly_activity(Ass_pos,data_h_rem');



%%
all_proj_rem = [];
for ii = size(wake_time_proj_rem, 1):-1:1
    
    %    nan_idx = time_proj_pos(ii,:) < prctile(time_proj_pos(ii,:), 10);
    %    plot3(ms_trk.time(~nan_idx)/1000, zscore(time_proj_pos(ii,~nan_idx)),ones(1, sum(~nan_idx))*100* ii, '.', 'color', c_ord(ii,:))
    %        plot3(tvec, zscore(time_proj_pos(ii,:)),ones(1, length(time_proj_pos(ii,:)))*100* ii, 'color', c_ord(ii,:))
    
    this_proj = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
    all_proj_rem(ii,:) =  this_proj;
    
end

%% split the REM assemblies based on location on the track.


close_idx = find(Ass_map_peak <= 50);
open_idx = find(Ass_map_peak >50);


close_peaks = sum(all_proj_rem(close_idx,:) > 1,2);
open_peaks = sum(all_proj_rem(open_idx,:) > 1,2);

% if sum(open_peaks == max(open_peaks)) ==1
[~, idx] = max(open_peaks);
Ass_1 = open_idx(idx(1));

[~, idx] = max(close_peaks);
Ass_2 = close_idx(idx(1));

% else % tie breaker
%
% end




%% collect the REM_react
rem_z.all = []; rem_z.close = []; rem_z.open = []; rem_z.isopen = [];
for ii =  size(wake_time_proj_rem,1):-1:1
    
    rem_z.all(ii,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
    
    if sum(ismember(close_idx, ii)) > 0
        rem_z.close(end+1,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
        rem_z.isopen(ii) = 0;
    elseif sum(ismember(open_idx, ii)) > 0
        rem_z.open(end+1,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
        rem_z.isopen(ii) = 1;
    end
end