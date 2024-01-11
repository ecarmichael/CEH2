function [rem_out, Ass_out, time_proj_pos_place, wake_time_proj_rem, behav_his] = MS_PCA_ICA_no_fig(fname)
% sandbox_PCA/ICA


if strcmp(computer, 'GLNXA64')
    
    codebase_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/home/williamslab/Documents/Github/CEH2';
    oasis_dir = '/home/williamslab/Documents/Github/OASIS_matlab';
    
    code_dir = '/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '/home/williamslab/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = '/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter';
    
    
    
else
    
    codebase_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
    ca_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
    
    code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';
    
    RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
    
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
plot_flag = 1;
% load('ms_trk.mat')
% load('behav_DLC.mat')
this_sess = fname;

load(this_sess)

behav = MS_align_data(behav,ms);

move_idx = behav.speed > 2.5;


behav_his = histcounts(behav.position(move_idx,1), 5:5:95);
%% check for place cell metrics
dir_parts = strsplit(this_sess, filesep);
parts = strsplit(dir_parts{end}, '_');
task = parts{2};
if contains(task, 'HATD6')
    task = 'HATDSwitch';
end
subject = parts{1};


if exist([subject '_' task '_PCs.mat'])
    load([subject '_' task '_PCs.mat'])
else
    %
    % if exist([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'], 'file')
    %     load([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'])
    % elseif exist([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'], 'file')
    %     load([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'])
    % else
    rem_out = NaN; Ass_out = NaN;  time_proj_pos_place = NaN; wake_time_proj_rem = NaN;
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
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk);

if ~isfield(ms_trk, 'deconv') % get the deconvolved trace if not already present.
    ms_trk_cut = MS_append_deconv(ms_trk_cut, 1);
end
% remove inactive cells

keep_idx = sum(ms_trk_cut.deconv, 1) >0;

remove_cell_id_decon = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv';
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk_cut);

%% grab the spatial tuning properties.


place = [];
for iC = length(PCs_properties.peak_loc):-1:1
    
    place.centroids(iC) = PCs_properties.peak_loc(iC);
    
    % is it a place cell?
    place.is(iC) = PCs_properties.isPC(iC);
    
    place.map(iC,:) = PCs_properties.tuning_curve_data(:,iC)';
    
    %     bin_vect = 0:2.5:100;
    %     p_idx = find(bin_vect == place.centroids(iC));
    %
    %     place.map(iC,:) = zeros(length(bin_vect),1); %mean(spatial_analysis.bin{iC,1}.PlaceField,1)/max(mean(spatial_analysis.bin{iC,1}.PlaceField,1)); % get the 1d place field and normalize.
    %     place.map(iC,p_idx) = 1;
end

place.centroids(remove_cell_id)= [];
place.centroids(remove_cell_id_decon)= [];

place.is(remove_cell_id)= [];
place.is(remove_cell_id_decon)= [];

place.map(remove_cell_id,:)= [];
place.map(remove_cell_id_decon,:)= [];

bin = 3; %80/size(place.map,2);
p_bins = 0:bin:100;
p_bins = p_bins(1:end)+bin/2;
% see if there are any anxiety cells

[~,p_sort] = sort(place.centroids);


if plot_flag
    figure(300);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    subplot(5,1,1)
    histogram(behav.position(:,1),p_bins, 'Normalization', 'probability')
    xlim([p_bins(1) p_bins(end)]); title('Occupancy');
    
    subplot(5,1,2)
    histogram(place.centroids(logical(place.is)),p_bins,'Normalization', 'probability', 'FaceColor', 'r')
    xlim([p_bins(1) p_bins(end)]); title('Centroids');
    
    
    subplot(5,1,3:5)
    this_map = place.map(p_sort,:);
    this_map(~place.is,:) = [];
    imagesc(p_bins, 1:length(place.centroids(logical(place.is))),  this_map./max(this_map,[],2))
    xlim([p_bins(1) p_bins(end)]);
    xlabel('Location on track (cm)')
    ylabel('Cell ID')
end
%% follow grosmark et al. method of deconv preprocessing
Csp = ms_trk_cut.deconv./ms_trk_cut.denoise;
Csp = Csp > 0.01;
ms_trk_cut.Csp = Csp;

% cfg_plot.Ca_type = 'RawTraces';
% cfg_plot.plot_type = '2d';
% MS_plot_ca(cfg_plot, ms_trk_rem)

% %% bin and convolve
% binsize = 0.1; % in seconds, so everything else should be seconds too
% gauss_window = 1./binsize; % 1 second window
% gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
% gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
% gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window
%
% gau_z = zscore(gau_sdf, [], 2);
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
    
    %     this_cell = ms.time(find(Csp(: ,ii) & move_idx))/1000;
    this_cell = ms.time(find(ms_trk_cut.Binary(: ,ii) & move_idx))/1000;
    
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
nShuff = 50;

Ass_shuff = NaN(1,nShuff);
for iS = 1:nShuff
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

Ass_out = (size(Ass_pos,2) - mean(Ass_shuff))/std(Ass_shuff);
fprintf('%.0f Positive Assemblies detected. Chance level is %.1f zscore = %.1fSD\n', size(Ass_pos,2), mean(Ass_shuff), Ass_out)
%% stem plot ensembles

Ass_map = cell(size(Ass_pos,2),1);
Ass_mean_map = [];
Ass_pcells = Ass_map;
c_ord = parula(size(Ass_pos,2)+2);
%


if plot_flag
    figure(301);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    cmap = parula(256);
    cmap(1,:) = 0;
end


for ii = size(Ass_pos,2):-1:1
    if plot_flag
        figure(301)
        subplot(4, ceil(size(Ass_pos,2)/4),ii)
        hold on
        stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
        view(90,90)
        
        stem(Ass_pos_cells{ii}, Ass_pos(Ass_pos_cells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
        title(['Assembly #' num2str(ii)])
        ylim([-.2 .4])
        
        
        %             figure(302)
        %                 subplot(4, ceil(size(Ass_pos,2)/4),ii)
        %             hold on
    end
    this_ass_map = []; these_place = zeros(size(Ass_pos_cells{ii}));
    fprintf('Assembly # %.0f had %.0f place cells\n', ii, sum(place.is(Ass_pos_cells{ii})))
    for jj = 1:length(Ass_pos_cells{ii})
        
        if place.is(Ass_pos_cells{ii}(jj))
            these_place(jj) = 1;
            place_int = interp1(p_bins,place.map(Ass_pos_cells{ii}(jj),:),  p_bins(1):1:p_bins(end));
            %             place_int = place.map(Ass_pos_cells{ii}(jj),:);
            this_ass_map = [this_ass_map ; place_int];
        end
    end
    
    if ~isempty(this_ass_map)
        %         if plot_flag
        %             imagesc(p_bins(1):1:p_bins(end), 1:length(sum(these_place)),  this_ass_map)
        %             set(gca,'ytick',1:size(this_ass_map,1), 'YTickLabel',  Ass_pos_cells{ii}(logical(these_place))');
        %             xlim([p_bins(1) p_bins(end)])
        %             ylim([.5 size(this_ass_map,1)+.5])
        %             colormap(cmap)
        %         end
        
        Ass_map{ii} = this_ass_map;
        Ass_pcells{ii} = Ass_pos_cells{ii}(logical(these_place));
        Ass_mean_map(ii,:) = mean(this_ass_map./max(this_ass_map,[],2),1);
        
        figure(302)
        subplot(4, ceil(size(Ass_pos,2)/4),ii)
        imagesc(p_bins(1):1:p_bins(end), 1:length(Ass_pcells{ii}),  Ass_map{ii}./max(Ass_map{ii},[],2))
        
    else
        
        Ass_map{ii} = NaN(size(p_bins(1):1:p_bins(end)));
        Ass_pcells{ii} = NaN;
        Ass_mean_map(ii,:) = NaN(size(p_bins(1):1:p_bins(end)));
    end
    
    %     if plot_flag
    %         title(['Assembly #' num2str(ii)])
    %         ms_t = MS_append_sharp_SFPs(ms_trk);
    %         MS_plot_all_SFPs(imgaussfilt(ms_t.SFPs_sharp(:,:, Ass_pos_cells{ii}),2))
    %     end
    
    
end
%%
nan_idx = isnan(max(Ass_mean_map, [], 2));
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



rm_idx =  find(cellfun(@length, Ass_pcells) < 2);
Ass_map_peak(rm_idx) = [];
Ass_idx(rm_idx) = [];
Ass_pos(:,rm_idx) = [];
Ass_mean_map(rm_idx,:) = [];
Ass_map(rm_idx) = [];
Ass_pcells(rm_idx) = [];
Ass_pos_cells_place(rm_idx) = [];
time_proj_pos_place(rm_idx,:) = [];



%% get the number of significant reactivations during wake for

c_ord = MS_linspecer(size(time_proj_pos_place,1)+4);

Wake_react = [];
Ass_rect_idx = [];
for ii = size(time_proj_pos_place,1):-1:1
    
    [this_rec, Ass_rect_idx{ii}] = findpeaks(time_proj_pos_place(ii,:), 'MinPeakHeight', 5, 'MinPeakDistance', 2/binsize);
    fprintf('Assembly #%d - %.0f sig reactivations (%0.2f/min)\n', ii, length(this_rec), length(this_rec)/((tvec(end)- tvec(1))/60))
    
end


if plot_flag
    figure(303);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    ax(1) = subplot(4,1,1);
    scatter(behav.time/1000, behav.position(:,1),ones(size(behav.time)), behav.speed)
    xlim([behav.time(1) behav.time(end)]/1000)
    colorbar('northoutside')
    
    ax(2) = subplot(4,1,2:4);
    hold on
    for  ii = size(time_proj_pos_place,1):-1:1
        plot(tvec, time_proj_pos_place(ii,:), 'color', c_ord(ii,:))
        
    end
    linkaxes(ax, 'x')
    
    %      figure(304);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %
    %      t_win = ceil(2/behav.dt);
    %
    %      for ii = length(Ass_rect_idx):-1:1
    %          subplot(4, ceil(size(Ass_rect_idx,2)/4),ii)
    %          hold on
    %
    %          this_idx = nearest_idx(tbin_centers(Ass_rect_idx{ii}), behav.time/1000);
    %
    %          for jj = length(this_idx):-1:1
    %         	plot((behav.time(this_idx(jj)-t_win:this_idx(jj)+t_win,1) - behav.time(this_idx(jj)-t_win))/1000, behav.position(this_idx(jj)-t_win:this_idx(jj)+t_win,1)', 'color', c_ord(ii,:))
    %
    %          end
    %
    %
    %      end
    
end

%% plot the remaining assembly maps;

if plot_flag
    cmap = parula(64);
    cmap(1,:) = 0;
    
    figure(302);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    figure(3021);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    for ii = size(Ass_pos,2):-1:1
        figure(302)
        subplot(4, ceil(size(Ass_pos,2)/4),ii)
        hold on
        stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
        view(90,90)
        
        stem(Ass_pcells{ii}, Ass_pos(Ass_pcells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
        title(['Assembly #' num2str(ii)])
        ylim([-0.1 0.4])
        
        
        fprintf('Assembly # %.0f had %.0f place cells\n', ii, sum(place.is(Ass_pos_cells_place{ii})))
        figure(3021)
        
        subplot(4, ceil(size(Ass_pos,2)/4),ii)
        %             imagesc(p_bins(1):1:p_bins(end), 1:length(Ass_pcells{ii})+1,  [Ass_map{ii}; nan(size(ass Ass_mean_map(ii,:)] )
        imagesc(p_bins(1):1:p_bins(end), 1:length(Ass_pcells{ii}),  Ass_map{ii}./max(Ass_map{ii}, [],2))
        
        set(gca,'ytick',1:length(Ass_pcells{ii}), 'YTickLabel',  Ass_pcells{ii}');
        xlim([p_bins(1) p_bins(end)])
        ylim([.5 size(Ass_map{ii},1)+.5])
        colormap(cmap)
        
    end
end



if plot_flag
    
    win = floor(2.5 * mode(diff(behav.time)));
    figure(305);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    clf
    maximize
    n = ceil(size(time_proj_pos_place,1)/3);
    m = 3;
    c_ord = MS_linspecer(size(time_proj_pos_place,1));
    
    for ii = 1: size(time_proj_pos_place,1)
        subplot(m,n,ii)
        hold on
        [p_val, p_idx] = findpeaks((time_proj_pos_place(ii,:)),'MinPeakHeight', 5 ,'MinPeakDistance', 2*floor(mode(diff(ms.time))));
        
        this_pos = [];
        for ip = 1:length(p_idx)
            this_idx = nearest_idx(tbin_centers(p_idx(ip)), behav.time/1000);
            if ((this_idx - win) >=0) && ((this_idx + win)<= length(behav.time))
                this_pos(ip,:) = behav.position(this_idx - win:this_idx+win,1);
                plot((-win:win)/mode(diff(behav.time)), this_pos(ip,:), 'color',[c_ord(ii,:) .5], 'linewidth',  2*(p_val(ip)./max(p_val)))
            end
        end
        
        plot((-win:win)/mode(diff(behav.time)), mean(this_pos), 'color',[c_ord(ii,:) 1], 'linewidth', 3)
        xlim([-win/mode(diff(behav.time)) win/mode(diff(behav.time))]);
        %         set(gca, 'color', 'k')
        title(['A#' num2str(ii)])
        plot(0, mean(this_pos(:,win)), 's','color','k', 'markersize', 20 )
        
        if ii == n+1
            xlabel({'time from assembly' ;  'onset (s)'})
        end
        if ii == 1 || ii== n+1
            ylabel('position on track (cm)')
        end
    end
end
%% summary figure


if plot_flag
    figure(3000);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    win = floor(2.5 * mode(diff(behav.time)));
    
    n = ceil(size(time_proj_pos_place,1)/3);
    m = 4;
    l = 5;
    c_ii = 0;
    f_n = 0;
    s_plot_max = reshape(1:(l*m), m, l)';
    
    for ii = 1: size(time_proj_pos_place,1)
        if c_ii+1 > l
            f_n = f_n+1;
            figure(3000+f_n);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            c_ii = 1;
        else
            c_ii = c_ii+1;
        end
        
        s_idx = s_plot_max(c_ii,1);
        
        subplot(l,m,s_idx)
        hold on
        stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
        view(90,90)
        
        stem(Ass_pcells{ii}, Ass_pos(Ass_pcells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
        title(['Assembly #' num2str(ii) ' (' num2str(length(Ass_pcells{ii})) ' place cells)'])
        ylim([-0.1 0.4])
        ylabel('cell ID')
        xlim([0 length(Ass_pos(:,ii))])
        
        
        subplot(l,m,s_idx+1)
        hold on
        [p_val, p_idx] = findpeaks((time_proj_pos_place(ii,:)),'MinPeakHeight', 10 ,'MinPeakDistance', 2*floor(mode(diff(ms.time))));
        
        this_pos = [];
        for ip = 1:length(p_idx)
            this_idx = nearest_idx(tbin_centers(p_idx(ip)), behav.time/1000);
            if ((this_idx - win) >=0) && ((this_idx + win)<= length(behav.time))
                this_pos(ip,:) = behav.position(this_idx - win:this_idx+win,1);
                plot((-win:win)/mode(diff(behav.time)), this_pos(ip,:), 'color',[c_ord(ii,:) .5], 'linewidth',  2*(p_val(ip)./max(p_val)))
            end
        end
        
        plot((-win:win)/mode(diff(behav.time)), median(this_pos), 'color',[c_ord(ii,:) 1], 'linewidth', 3)
        xlim([-win/mode(diff(behav.time)) win/mode(diff(behav.time))]);
        %         set(gca, 'color', 'k')
        plot(0, median(this_pos(:,win)), 's','color','k', 'markersize', 20 )
        xlabel('time from ReAct (s)')
        ylabel('position on track (cm)')
        
        
        subplot(l,m,s_idx+2)
        cla
        hold on
        [N, edges] = histcounts(this_pos, p_bins(1):3:p_bins(end));
        
        area(p_bins(1):1:p_bins(end),Ass_mean_map(ii,:)./max( Ass_mean_map(ii,:)), 'facecolor',c_ord(ii,:), 'EdgeAlpha', 0 )
        plot(edges(1:end-1)+mode(diff(edges))/2, N./max(N),'-', 'color', [0.3 .3 .3], 'linewidth', 1)
        set(gca, 'YTick',[0 1], 'yticklabel', {'0' 'max'}, 'xdir', 'reverse')
        xlim([p_bins(1) p_bins(end)])
        ylabel('Mean place map')
        view(90, 90)
        legend({'Mean place field', 'ReAct location'}, 'Location', 'northeast', 'Box', 'off')
        
        
        
        subplot(l,m,s_idx+3)
        %             imagesc(p_bins(1):1:p_bins(end), 1:length(Ass_pcells{ii})+1,  [Ass_map{ii}; nan(size(ass Ass_mean_map(ii,:)] )
        imagesc(1:length(Ass_pcells{ii}),p_bins(1):1:p_bins(end),  (Ass_map{ii}./max(Ass_map{ii}, [],2))')
        
        set(gca,'xtick',1:length(Ass_pcells{ii}), 'xTickLabel',  Ass_pcells{ii}','ydir', 'normal');
        ylim([p_bins(1) p_bins(end)])
        xlim([.5 size(Ass_map{ii},1)+.5])
        colormap(cmap)
        xlabel('place cell ID')
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grab the REM data

all_binary_post_REM(:, remove_cell_id) = [];
all_binary_post_REM(:, remove_cell_id_decon) = [];

% all_deconv_post_REM =


all_binary_pre_REM(:, remove_cell_id) = [];
all_binary_pre_REM(:, remove_cell_id_decon) = [];


%% deconvolve
if exist('all_detrendRaw_post_REM')
    all_detrendRaw_post_REM(:, remove_cell_id) = [];
    all_detrendRaw_post_REM(:, remove_cell_id_decon) = [];
    
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
    
end
%% Get the Post REM activity
this_time =  0:1/mode(diff(ms.time)):(length(all_binary_post_REM)/mode(diff(ms.time)));
this_time = this_time(1:end-1);

binsize = .5;
tbin_edges_rem = 0:binsize:(length(all_binary_post_REM)/mode(diff(ms.time))); % vector of time bin edges (for histogram)
tbin_centers_rem = tbin_edges_rem(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_h_rem = [];
for ii = size(all_binary_post_REM,2):-1:1
    
    if exist('all_detrendRaw_post_REM', 'var')
        this_cell = this_time(find(Csp_rem(: ,ii)));
    else
        this_cell = this_time(find(all_binary_post_REM(:,ii)));
    end
    spk_count = histc(this_cell,tbin_edges_rem); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    data_h_rem(:,ii) = spk_count;
end
tvec_rem = tbin_centers_rem;

%% Get the Pre REM activity
this_time =  0:1/mode(diff(ms.time)):(length(all_binary_pre_REM)/mode(diff(ms.time)));
this_time = this_time(1:end-1);

binsize = .5;
tbin_edges_rem = 0:binsize:(length(all_binary_pre_REM)/mode(diff(ms.time))); % vector of time bin edges (for histogram)
tbin_centers_rem = tbin_edges_rem(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_h_rem_pre = [];
for ii = size(all_binary_pre_REM,2):-1:1
    
    if exist('all_detrendRaw_post_REM', 'var')
        this_cell = this_time(find(Csp_rem(: ,ii)));
    else
        this_cell = this_time(find(all_binary_pre_REM(:,ii)));
    end
    spk_count = histc(this_cell,tbin_edges_rem); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    data_h_rem_pre(:,ii) = spk_count;
end
tvec_rem_pre = tbin_centers_rem;

%% try the assembly code
rng(123, 'twister')
Ass_Temp_rem = assembly_patterns(data_h_rem');

rng(123, 'twister')
%using wake assemblies.
wake_time_proj_rem = assembly_activity(Ass_pos,data_h_rem');

% PRE using wake assemblies.
rng(123, 'twister')
wake_time_proj_rem_pre = assembly_activity(Ass_pos,data_h_rem_pre');

rng(123, 'twister')
all_time_proj_rem = assembly_activity(Ass_Temp,data_h_rem');

sig_REM_react = [];
for ii = size(wake_time_proj_rem,1):-1:1
    [~, p_idx] = findpeaks((wake_time_proj_rem(ii,:)),'MinPeakHeight', 10 ,'MinPeakDistance', 2/(mode(diff(tvec_rem))));
    
    if ~isempty(p_idx)
        sig_REM_react(ii) = length(p_idx);
    else
        sig_REM_react(ii) = NaN;
    end
    
end


for ii = 1:size(Ass_pos,2)
        
        
        % pre assembly react mean
        [~, pre_idx] = findpeaks(wake_time_proj_rem_pre(ii,:),'MinPeakHeight', 20 ,'MinPeakDistance', 2/(mode(diff(tvec_rem_pre))));
        [~, post_idx] = findpeaks(wake_time_proj_rem(ii,:),'MinPeakHeight', 20 ,'MinPeakDistance', 2/(mode(diff(tvec_rem))));

        fprintf('Assembly # %f had PRE %f  & POST %f sig reactivations (threst > 20)\n', ii, length(pre_idx), length(post_idx))
end

%% get the REM react shuffle.

rng(123,'twister')
nShuff = 500;
shuff_mat = []; 
shuff_mat_pre = []; 

for iS = 1:nShuff
    tic
    shuff_data = NaN(size(data_h_rem_pre));
    for ic = 1:size(data_h_rem_pre,2)
        shuff_data(:,ic) = circshift(data_h_rem_pre(:,ic), floor(MS_randn_range(1,1,1,size(data_h_rem_pre,1))));
    end
    
    wake_time_proj_rem_s = assembly_activity(Ass_pos,shuff_data');
    
    shuff_mat_pre = [shuff_mat_pre; wake_time_proj_rem_s];
    
    
   shuff_data = NaN(size(data_h_rem));
    for ic = 1:size(data_h_rem,2)
        shuff_data(:,ic) = circshift(data_h_rem(:,ic), floor(MS_randn_range(1,1,1,size(data_h_rem,1))));
    end
    
    wake_time_proj_rem_s = assembly_activity(Ass_pos,shuff_data');
    
    shuff_mat = [shuff_mat; wake_time_proj_rem_s];

end

R_threshold = prctile(shuff_mat(shuff_mat >0), 99, 'all'); 

%% loop over assemblies to see what the range of reactivations would be. 
Ass_p_val_pre = []; 
ReAct_rate_pre = []; 
ReAct_rate_p_pre = []; 


Ass_p_val_post = []; 
ReAct_rate_post = []; 
ReAct_rate_p_post = []; 


for ii = size(Ass_pos,2):-1:1    
    Ass_p_val_pre(ii) = sum(sum(shuff_mat_pre(1:1000,:)>R_threshold) >sum(wake_time_proj_rem_pre(ii,:) > R_threshold))/ numel(shuff_mat_pre(1:1000,:));
    
    ReAct_rate_pre(ii) = sum(wake_time_proj_rem_pre(ii,:) > R_threshold) / ((tvec_rem_pre(end) - tvec_rem(1))/60); 

    Shuff_rate_pre = sum(shuff_mat_pre(1:1000,:) > R_threshold,2)./ ((tvec_rem_pre(end) - tvec_rem_pre(1))/60); 
    
    ReAct_rate_p_pre(ii) = sum(Shuff_rate_pre > ReAct_rate_pre(ii)) / length(Shuff_rate_pre); 
    
    
    
    Ass_p_val_post(ii) = sum(sum(shuff_mat(1:1000,:)>R_threshold) >sum(wake_time_proj_rem(ii,:) > R_threshold))/ numel(shuff_mat(1:1000,:));
    
    ReAct_rate_post(ii) = sum(wake_time_proj_rem(ii,:) > R_threshold) / ((tvec_rem(end) - tvec_rem(1))/60); 

    Shuff_rate_post = sum(shuff_mat(1:1000,:) > R_threshold,2)./ ((tvec_rem(end) - tvec_rem(1))/60); 
    
    ReAct_rate_p_post(ii) = sum(Shuff_rate_post > ReAct_rate_post(ii)) / length(Shuff_rate_post); 

end



% shuff_react = [];
% for ii = size(wake_time_proj_rem_s,1):-1:1
%     [~, p_idx] = findpeaks((wake_time_proj_rem_s(ii,:)),'MinPeakHeight', 20 ,'MinPeakDistance', 2/(mode(diff(tvec_rem))));
%     if ~isempty(p_idx)
%         shuff_react(ii) = length(p_idx);
%         
%     else
%         shuff_react(ii) = NaN;
%     end
% end
% shuff_sig_react(iS) = sum(~isnan(shuff_react));
% 
% fprintf('Shuff # %.0f found %.0f assemblies and took %2.2f seconds\n', iS, shuff_react(ii), toc)



if plot_flag
    figure(999); clf; 
    histogram(shuff_mat); 

end


rem_out.shuff_time_prog_rem_z = shuff_time_prog_rem_z;

%%
if plot_flag
    A_idx = 1:size(Ass_pos,2);
    %     A_idx = [3, 4, 9 12 13];
    figure(310);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    ar(1)= subplot(6,1,1);
    
    imagesc(tvec_rem_pre,1:size(data_h_rem_pre,2),   data_h_rem_pre')
    xlim([tvec_rem_pre(1) tvec_rem_pre(end)])
    title([subject ' ' task])
    set(gca, 'xtick', [])
    
    
    ar(2) = subplot(6,1,2);
    plot(tvec_rem_pre, wake_time_proj_rem_pre(A_idx,:))
    ylim([0 50])
    xlim([tvec_rem_pre(1) tvec_rem_pre(end)])
    ylabel('Pre REM Reactivation')
    legend(num2str(A_idx'),'location',  'northeast', 'box', 'off')
    
    ap(1) = subplot(6,1,3);
    imagesc(tvec_rem,1:size(data_h_rem,2),   data_h_rem')
    xlim([tvec_rem(1) tvec_rem(end)])
        set(gca, 'xtick', [])

    ap(2) = subplot(6,1,4);
    plot(tvec_rem, wake_time_proj_rem(A_idx,:))
    ylim([0 50])
    xlim([tvec_rem(1) tvec_rem(end)])
    ylabel('Post REM Reactivation')
    
    as(1) = subplot(6,1,5);
    imagesc(tvec_rem,1:size(shuff_data,2),   shuff_data')
    xlim([tvec_rem(1) tvec_rem(end)])
        set(gca, 'xtick', [])

    as(2) = subplot(6,1,6);
    plot(tvec_rem, wake_time_proj_rem_s(A_idx,:))
    ylim([0 50])
    xlim([tvec_rem(1) tvec_rem(end)])
    ylabel('Shuff postReactivation')
    xlabel('time (s)')
    
    
    linkaxes(ar, 'x');
    linkaxes(ap, 'x');
    linkaxes(as, 'x');
    
end
% %%
% all_proj_rem = [];
% for ii = size(wake_time_proj_rem, 1):-1:1
%
%     %    nan_idx = time_proj_pos(ii,:) < prctile(time_proj_pos(ii,:), 10);
%     %    plot3(ms_trk.time(~nan_idx)/1000, zscore(time_proj_pos(ii,~nan_idx)),ones(1, sum(~nan_idx))*100* ii, '.', 'color', c_ord(ii,:))
%     %        plot3(tvec, zscore(time_proj_pos(ii,:)),ones(1, length(time_proj_pos(ii,:)))*100* ii, 'color', c_ord(ii,:))
%
%     this_proj = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
%     all_proj_rem(ii,:) =  this_proj;
%
% end
%% react triggered spike patterns
if plot_flag
    figure(311);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    m = 4; n = 6;  % for subplots
    c = repelem(1:2:(n*m)/2, 2);
    %     win = floor(.2 * mode(diff(behav.time)));
    win = 1;
    for ii = 1:size(Ass_pos,2)
        
        
        % pre assembly react mean
        [~, pre_idx] = findpeaks(wake_time_proj_rem_pre(ii,:),'MinPeakHeight', 20 ,'MinPeakDistance', 2/(mode(diff(tvec_rem_pre))));
        
        subplot(m, n, c(ii))
        if ~isempty(pre_idx)
            this_mat = nan(length(-win:win), size(data_h_rem_pre,2), length(pre_idx));
            for jj = length(pre_idx):-1:1
                
                if min(pre_idx(jj)-win:pre_idx(jj)+win) <1
                    z_idx = find(pre_idx(jj)-win:pre_idx(jj)+win ==1);
                    l = length(-win:win);
                    this_mat(z_idx:l,:,jj) = data_h_rem_pre(pre_idx(jj)-(win-z_idx+1):pre_idx(jj)+win,:);
                    
                elseif max(pre_idx(jj)-win:pre_idx(jj)+win) > size(data_h_rem_pre,1)
                    z_idx = find(pre_idx(jj)-win:pre_idx(jj)+win ==size(data_h_rem_pre,1));
                    l = length(-win:win);
                    this_mat(z_idx:l,:,jj) = data_h_rem_pre(pre_idx(jj)-win:pre_idx(jj)+win,:);
                else
                    this_mat(:,:,jj) = data_h_rem_pre(pre_idx(jj)-win:pre_idx(jj)+win,:);
                end
            end
            
            imagesc((-win:win)/mode(diff(behav.time)), 1:size(this_mat,2), nanmean(this_mat, 3)');
            set(gca, 'ydir', 'normal')
            title(['Assembly ' num2str(ii) ' Pre'])
        end
        
        
        
        % post assembly react mean
        [~, post_idx] = findpeaks(wake_time_proj_rem(ii,:),'MinPeakHeight', 20 ,'MinPeakDistance', 2/(mode(diff(tvec_rem))));
        
        subplot(m, n, c(ii)+1)
        if ~isempty(post_idx)
            this_mat = nan(length(-win:win), size(data_h_rem,2), length(post_idx));
            for jj = length(post_idx):-1:1
                
                if min(post_idx(jj)-win:post_idx(jj)+win) <1
                    z_idx = find(post_idx(jj)-win:post_idx(jj)+win ==1);
                    l = length(-win:win);
                    this_mat(z_idx:l,:,jj) = data_h_rem(post_idx(jj)-(win-z_idx+1):post_idx(jj)+win,:);
                    
                elseif max(post_idx(jj)-win:post_idx(jj)+win) > size(data_h_rem,1)
                    z_idx = find(post_idx(jj)-win:post_idx(jj)+win ==size(data_h_rem,1));
                    this_mat(end-z_idx+1:end,:,jj) = data_h_rem(post_idx(jj)-win:end,:);
                else
                    this_mat(:,:,jj) = data_h_rem(post_idx(jj)-win:post_idx(jj)+win,:);
                end
            end
            imagesc((-win:win)/mode(diff(behav.time)), 1:size(this_mat,2), nanmean(this_mat, 3)');
            set(gca, 'ydir', 'normal')
            title(['Assembly ' num2str(ii) ' Post'])
        end
        
    end
end


%% split the REM assemblies based on location on the track.
if ~isempty(Ass_map_peak)
    
    
    if strcmpi(task, 'HATDS')
        close_idx = find(Ass_map_peak < 40);
        mid_idx = find((40 <= Ass_map_peak) & (Ass_map_peak<=60));
        open_idx = find(Ass_map_peak >60);
    else
        close_idx = find(Ass_map_peak > 60);
        mid_idx = find((40 <= Ass_map_peak) & (Ass_map_peak<=60));
        
        open_idx = find(Ass_map_peak <40);
    end
    
    
    % close_peaks = sum(all_proj_rem(close_idx,:) > .4,2);
    % open_peaks = sum(all_proj_rem(open_idx,:) > .4,2);
    
    % if sum(open_peaks == max(open_peaks)) ==1
    % [~, idx] = max(open_peaks);
    % Ass_1 = open_idx(idx(1));
    %
    % [~, idx] = max(close_peaks);
    % Ass_2 = close_idx(idx(1));
else
    rem_out = NaN;
    return
end

% else % tie breaker
%
% end




%% collect the REM_react
rem_out.all = []; rem_out.close = []; rem_out.open = []; rem_out.mid = []; rem_out.isopen = [];rem_out.ismid = []; rem_out.close_z = []; rem_out.open_z = []; rem_out.mid_z = [];
rem_out.all_time_prog = all_time_proj_rem;
rem_out.all_time_prog_z = all_time_prog_rem_z;

for ii =  size(wake_time_proj_rem,1):-1:1
    
    rem_out.all(ii,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
    
    if sum(ismember(close_idx, ii)) > 0
        rem_out.close(end+1,:) = wake_time_proj_rem(ii,:);
        rem_out.close_z(end+1,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
        
        rem_out.isopen(ii) = 0;
        rem_out.ismid(ii) = 0;
        
    elseif sum(ismember(open_idx, ii)) > 0
        rem_out.open_z(end+1,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
        rem_out.open(end+1,:) = wake_time_proj_rem(ii,:);
        rem_out.isopen(ii) = 1;
        rem_out.ismid(ii) = 0;
        
    elseif sum(ismember(mid_idx, ii)) > 0
        rem_out.mid_z(end+1,:) = wake_time_proj_rem(ii,:);
        rem_out.mid(end+1,:) = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:));
        
        rem_out.ismid(ii) = 1;
        rem_out.isopen(ii) = 0;
        
    end
end