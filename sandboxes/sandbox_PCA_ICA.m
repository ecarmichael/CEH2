function [rem_z, Ass_z] = sandbox_PCA_ICA(fname)
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
    rem_z = NaN; Ass_z = NaN; 
    return
end
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
figure(301);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
figure(302);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

Ass_map = cell(size(Ass_pos,2),1); 
Ass_mean_map = []; 
Ass_pcells = Ass_map;
c_ord = parula(size(Ass_pos,2)+2);

cmap = parula(256);
cmap(1,:) = 0; 
for ii = size(Ass_pos,2):-1:1
    figure(301)
    subplot(4, ceil(size(Ass_pos,2)/4),ii)
    hold on
    stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
    view(90,90)
    
    stem(Ass_pos_cells{ii}, Ass_pos(Ass_pos_cells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
        title(['Assembly #' num2str(ii)])

    
    figure(302)
        subplot(4, ceil(size(Ass_pos,2)/4),ii)
    hold on
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
        imagesc(p_bins(1):1:p_bins(end), 1:length(sum(these_place)),  this_ass_map)
        set(gca,'ytick',1:size(this_ass_map,1), 'YTickLabel',  Ass_pos_cells{ii}(logical(these_place))');
        xlim([p_bins(1) p_bins(end)])
        ylim([.5 size(this_ass_map,1)+.5])
        colormap(cmap)
        Ass_map{ii} = this_ass_map;
        Ass_pcells{ii} = Ass_pos_cells{ii}(logical(these_place));
        Ass_mean_map(ii,:) = mean(this_ass_map,1);
    else
        
      Ass_map{ii} = NaN(size(p_bins(1):1:p_bins(end))); 
      Ass_pcells{ii} = NaN; 
      Ass_mean_map(ii,:) = NaN(size(p_bins(1):1:p_bins(end))); 
    end
    title(['Assembly #' num2str(ii)])
    


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
figure(1001); clf; 
[~, idx] = max(Ass_mean_map, [], 2); 
p_bins_int = p_bins(1):1:p_bins(end); 

Ass_map_peak = p_bins_int(idx); 
[Ass_map_peak, s_idx] = sort(Ass_map_peak); 
subplot(1,3,1)
imagesc(p_bins_int, 1:size(Ass_mean_map,1), Ass_mean_map);
set(gca,'ytick',1:size(Ass_mean_map,1)-.5,  'YTickLabel', num2str(Ass_idx(:)))

Ass_idx = Ass_idx(s_idx); 
Ass_pos = Ass_pos(:,s_idx); 
Ass_mean_map = Ass_mean_map(s_idx,:); 
Ass_map = Ass_map(s_idx); 
Ass_pcells = Ass_pcells(s_idx); 
Ass_pos_cells_place = Ass_pos_cells_place(s_idx); 
time_proj_pos_place = time_proj_pos_place(s_idx,:); 

subplot(1,3,2)
imagesc(p_bins_int, 1:size(Ass_mean_map,1), Ass_mean_map);
set(gca,'ytick',1:size(Ass_mean_map,1)-.5,  'YTickLabel', num2str(Ass_idx(:)))

rm_idx =  find(cellfun(@length, Ass_pcells) < 3); 
Ass_map_peak(rm_idx) = [];
Ass_idx(rm_idx) = []; 
Ass_pos(:,rm_idx) = []; 
Ass_mean_map(rm_idx,:) = []; 
Ass_map(rm_idx) = [];
Ass_pcells(rm_idx) = []; 
Ass_pos_cells_place(rm_idx) = [];  
time_proj_pos_place(rm_idx,:) = []; 
subplot(1,3,3)
imagesc(p_bins_int, 1:size(Ass_mean_map,1), Ass_mean_map);
set(gca,'ytick',1:size(Ass_mean_map,1)-.5,  'YTickLabel', num2str(Ass_idx(:)))
%% updated the plots
figure(303);clf; hold on
figure(304);clf; hold on

c_ord = parula(size(Ass_pos,2)+2);


for ii = size(Ass_pos,2):-1:1
    figure(303)
    subplot(4, ceil(size(Ass_pos,2)/4),ii)
    hold on
    stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
    view(90,90)
    
    stem(Ass_pos_cells_place{ii}, Ass_pos(Ass_pos_cells_place{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
        title(['Assembly #' num2str(ii)])

    
    figure(304)
        subplot(4, ceil(size(Ass_pos,2)/4),ii)
    hold on
    
    imagesc(p_bins(1):1:p_bins(end), 1:size(Ass_map{ii},1), Ass_map{ii})
            xlim([p_bins(1) p_bins(end)])
        ylim([.5 size(Ass_map{ii},1)+.5])
                colormap(cmap)

        title(['Assembly#' num2str(Ass_idx(ii))])
end




%% get the number of significant reactivations during wake for 

%% get the mean pop activity as per SCE
data_in = Csp; 

df = mode(diff(ms.time)); 
frame_n = floor(df/2);
if mod(frame_n, 2)==0; frame_n = frame_n+1; end
frame_n_200 = floor(df/5); 

shuff = 100;
all_shuff = [];
tic
for iS = shuff:-1:1
    shuff_data = []; 
    for ii = size(data_in, 2):-1:1
        shuff_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
    end % end cells

    all_shuff(iS, :) = movmean(nansum(shuff_data,1), frame_n_200); 
end % end shuff
toc

shuff_mean = mean(all_shuff,'all');
shuff_sd = std(all_shuff, [],'all');

thresh = shuff_mean + 3*shuff_sd; 
% find time points in real data that exceed threshold

pop_act = movmean(nansum(data_in,2), frame_n_200);

% exlude events that are too close. use findpeaks
% figure(101); 
[peak_act, SCE_idx] = findpeaks(pop_act, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', df);
%% plot the output

figure(202)
clf
% maximize
ax(1) = subplot(6, 1, 1);
hold on
plot(ms_trk.time/1000, behav.position(:,1));
% plot(ms_trk.time/1000, behav.position(:,2));
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylabel('position on track (cm)')
set(gca, 'XTick', []);

ax(2) = subplot(6,1,2:3);
imagesc(ms.time/1000, 1:size(gau_z,2), gau_z')
% hold on
% for ii = size(gau_sdf, 2):-1:1
%     s_t = ms_trk.time(Csp(:,ii) >0)/1000;
%     if ~isempty(s_t)
%        plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', 'k', 'linewidth', 4)
%
%     end
% %     plot(ms_trk.time, gau_sdf(:,ii)+ii)
% end
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylim([0 size(Csp, 2)])
set(gca, 'XTick', [], 'YDir', 'normal')
ylabel('cell ID')

c = colorbar('Location', 'east');
c.AxisLocation = 'out';
pos = c.Position;
c.Position = [pos(1)+.025 pos(2) pos(3) pos(4)];
c.Label.String = 'z-score activity';
c.Label.FontSize = 20;

ax(3) = subplot(6,1,4:5);
cla
hold on
c_ord = cool(size(time_proj_pos_place,1)+4);
for ii = 1:size(time_proj_pos_place,1)
    %    nan_idx = time_proj_pos(ii,:) < prctile(time_proj_pos(ii,:), 10);
    %    plot3(ms_trk.time(~nan_idx)/1000, zscore(time_proj_pos(ii,~nan_idx)),ones(1, sum(~nan_idx))*100* ii, '.', 'color', c_ord(ii,:))
    %        plot3(tvec, zscore(time_proj_pos(ii,:)),ones(1, length(time_proj_pos(ii,:)))*100* ii, 'color', c_ord(ii,:))
    
%     if 
    
    plot(tvec, zscore(time_proj_pos_place(ii,:))+ii*10, 'color', c_ord(ii+4,:))
    
    tick_val(ii) = mode(zscore(time_proj_pos_place(ii,:))+ii*10);
    tick_label{ii} = num2str(ii);
end
% view(0,  15)
set(gca, 'YTick', tick_val, 'YTickLabel', tick_label)

ylabel('assembly ID')
xlabel('time (s)')
ylim([min(zscore(time_proj_pos_place(1,:))+10), max(zscore(time_proj_pos_place(1,:))+10*size(time_proj_pos_place,1))])
xlim([tvec(1) tvec(end)])
%     xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
%     zlim([0 size(time_proj_pos,1)*100])
%     set(gca, 'color', 'k')

% plot the SCE
ax(4) = subplot(6,1,6);
cla
hold on
plot(ms.time/1000, pop_act, 'color', 'k', 'linewidth', 2)
plot(ms.time/1000, nanmean(all_shuff) - shuff_sd*3,'--',  'color', [.8 .8 .8 .3], 'linewidth',.5)
plot(ms.time/1000, nanmean(all_shuff) + shuff_sd*3,'--', 'color', [.8 .8 .8 .3], 'linewidth',.5)
plot(ms.time/1000, nanmean(all_shuff), 'color', [.8 .8 .8 .3], 'linewidth', 1)

ylabel('pop activity')


linkaxes(ax, 'x')
xlim([ms.time(1)/1000 ms.time(end)/1000])


linkaxes(ax, 'x')




%% get the assembly triggered position average.
win = floor(2.5 * mode(diff(behav.time)));
figure(5)
clf
n = ceil(size(time_proj_pos_place,1)/2);
m = 2;

for ii = 1: size(time_proj_pos_place,1)
    subplot(m,n,ii)
    hold on
    [~, p_idx] = findpeaks(zscore(time_proj_pos_place(ii,:)),'MinPeakHeight', 1.96 ,'MinPeakDistance', 2*floor(mode(diff(ms.time))));
    
    this_pos = [];
    for ip = 1:length(p_idx)
        this_idx = nearest_idx(tbin_centers(p_idx(ip)), behav.time/1000);
        if ((this_idx - win) >=0) && ((this_idx + win)<= length(behav.time))
            this_pos(ip,:) = behav.position(this_idx - win:this_idx+win,1);
            plot((-win:win)/mode(diff(behav.time)), this_pos(ip,:), 'color',[c_ord(ii+4,:) .5])
        end
    end
    
    plot((-win:win)/mode(diff(behav.time)), mean(this_pos), 'color',[c_ord(ii+4,:) 1], 'linewidth', 3)
    xlim([-win/mode(diff(behav.time)) win/mode(diff(behav.time))]);
    %         set(gca, 'color', 'k')
    title(['Assembly #' num2str(ii)])
    plot(0, mean(this_pos(:,win)), 's','color',[c_ord(ii,:) 1], 'markersize', 20 )
    
    if ii == n+1
        xlabel({'time from assembly' ;  'onset (s)'})
    end
    if ii == 1 || ii== n+1
        ylabel('position on track (cm)')
    end
end


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

load('all_detrendRaw_post_REM.mat')
load('all_RawTraces_post_REM.mat')
load('all_binary_post_REM.mat')


%remove cells that were excluded in the awake set.

all_detrendRaw_post_REM(:, remove_cell_id) = [];
all_detrendRaw_post_REM(:, remove_cell_id_decon) = [];

all_RawTraces_post_REM(:, remove_cell_id) = [];
all_RawTraces_post_REM(:, remove_cell_id_decon) = [];

all_binary_post_REM(:, remove_cell_id) = [];
all_binary_post_REM(:, remove_cell_id_decon) = [];
% load the decoding as well

cd(decode_dir)

s_list = dir('pv*');

sessions = [];
if strcmpi(sess, 'HATSwitch')

sess = 'HATDSwitch';
end
keep_idx = zeros(1, length(s_list)); 
for ii = length(s_list):-1:1
if ~isempty(strfind(lower(s_list(ii).name), lower(sess))) && ~isempty(strfind(lower(s_list(ii).name), lower(sub)))
 keep_idx(ii) = 1; 
end
end

load(s_list(find(keep_idx)).name)

decode_prob = decoding.REM_decoded_probabilities; 
decode_pos = decoding.REM_decoded_position; 
decode_bin = decoding.bin_centers_vector; 
clear decoding; 

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


%% bin and convolve
binsize = 0.1; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf_rem = conv2(Csp_rem,gk,'same'); % convolve with gaussian window

gau_z_rem = zscore(gau_sdf_rem, [], 2);

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

%% get the mean pop activity as per SCE
data_in = all_binary_post_REM; 

df = mode(diff(rem_time*1000)); 
frame_n = floor(df/2);
if mod(frame_n, 2)==0; frame_n = frame_n+1; end
frame_n_200 = floor(df/5); 

shuff = 100;
all_shuff_rem = [];
tic
for iS = shuff:-1:1
    shuff_data = []; 
    for ii = size(data_in, 2):-1:1
        shuff_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
    end % end cells

    all_shuff_rem(iS, :) = movmean(nansum(shuff_data,1), frame_n_200); 
end % end shuff
toc

shuff_mean_rem = mean(all_shuff_rem,'all');
shuff_sd_rem = std(all_shuff_rem, [],'all');

thresh = shuff_mean + 3*shuff_sd_rem; 
% find time points in real data that exceed threshold

pop_act_rem = movmean(nansum(data_in,2), frame_n_200);

% exlude events that are too close. use findpeaks
% figure(101); 
[peak_act_rem, SCE_idx_rem] = findpeaks(pop_act_rem, 1, 'MinPeakHeight', thresh, 'MinPeakDistance', df);

%% try the assembly code
rng(123, 'twister')
Ass_Temp_rem = assembly_patterns(data_h_rem');

rng(123, 'twister')
%using wake assemblies.
wake_time_proj_rem = assembly_activity(Ass_pos,data_h_rem');

%% plot the wake Assemblies in REM

figure(5); clf;
ax(1) = subplot(6,1,1);
imagesc(rem_time, decode_bin,   decode_prob)
ylabel('decoded position')
set(gca, 'xtick', [], 'YDir', 'normal')
c = colorbar('Location', 'east');
c.AxisLocation = 'out';
pos = c.Position;
c.Position = [pos(1)+.025 pos(2) pos(3) pos(4)];
c.Label.String = {'decoded'; 'probability'};
c.Label.FontSize = 14;
caxis([0 .3])
title('Reactivation of wake assemblies in REM')


ax(2) = subplot(6,1,2:3);
cla
hold on
% offset = 10; 
% p_ord = parula(15); 
% for ii = 1:10
% %    plot(rem_time, zscore(all_detrendRaw_post_REM(:,ii))'+ii*offset); 
%    plot(rem_time, zscore(all_denoise(:,ii))'+ii*offset, 'color', p_ord(ii,:)); 
%    plot(rem_time, Csp_rem(:,ii)+(ii*offset) - offset/4,'color', [0 0 0])
%       plot(rem_time, all_binary_post_REM(:,ii)+(ii*offset) - offset/3,'color', [.5 .5 .5])
% 
%       plot(rem_time, zscore(all_RawTraces_post_REM(:,ii))'+ii*offset,'color', [p_ord(ii,:) .4]); 
% 
% 
% end
imagesc(tvec_rem, 1:size(all_binary_post_REM,2), all_binary_post_REM')
ylim([0.5 size(all_binary_post_REM,2)+.5])
c = colorbar('Location', 'east');
c.AxisLocation = 'out';
pos = c.Position;
c.Position = [pos(1)+.025 pos(2) pos(3) pos(4)];
c.Label.String = 'binarized activity';
c.Label.FontSize = 14;
% c.Ticks = 
ylabel('cell ID')
caxis([0 1])
set(gca, 'xtick', [], 'YDir', 'normal')



ax(3)= subplot(6,1,4:5);
cla
hold on
c_ord = cool(size(wake_time_proj_rem,1)+4);

for ii = 1:size(wake_time_proj_rem, 1)
    
    %    nan_idx = time_proj_pos(ii,:) < prctile(time_proj_pos(ii,:), 10);
    %    plot3(ms_trk.time(~nan_idx)/1000, zscore(time_proj_pos(ii,~nan_idx)),ones(1, sum(~nan_idx))*100* ii, '.', 'color', c_ord(ii,:))
    %        plot3(tvec, zscore(time_proj_pos(ii,:)),ones(1, length(time_proj_pos(ii,:)))*100* ii, 'color', c_ord(ii,:))
    
    this_proj = (wake_time_proj_rem(ii,:) - mean(time_proj_pos_place(ii,:)))/std(time_proj_pos_place(ii,:)); 
    all_proj_rem(ii,:) =  this_proj;
    plot(tvec_rem, this_proj+ii*1, 'color', c_ord(ii+4,:))
    
    tick_val(ii) = mode(this_proj)+ii;
    tick_label{ii} = num2str(ii);
end
% view(0,  15)
set(gca, 'YTick', tick_val(1:2:end), 'YTickLabel', tick_label(1:2:end))

ylabel('assembly ID')
xlabel('time (s)')
% ylim([min(zscore(wake_time_proj_rem(1,:))+10), max(zscore(wake_time_proj_rem(1,:))+10*size(wake_time_proj_rem,1))])
%     xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
%     zlim([0 size(time_proj_pos,1)*100])
%     set(gca, 'color', 'k')



% plot the SCE
ax(4) = subplot(6,1,6);
cla
hold on
plot(rem_time, pop_act_rem, 'color', 'k', 'linewidth', 2)
plot(rem_time, nanmean(all_shuff_rem) - shuff_sd_rem*3,'--',  'color', [.8 .8 .8 .3], 'linewidth',.5)
plot(rem_time, nanmean(all_shuff_rem) + shuff_sd_rem*3,'--', 'color', [.8 .8 .8 .3], 'linewidth',.5)
plot(rem_time, nanmean(all_shuff_rem), 'color', [.8 .8 .8 .3], 'linewidth', 1)

ylabel('pop activity')


linkaxes(ax, 'x')
xlim([rem_time(1) rem_time(end)])

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



%% break out figures for specific assemblies. 



% Ass_1 = 8; 
% Ass_2 = 1; 

Ass_1_ids = Ass_pos_cells_place{Ass_1}; 
Ass_2_ids = Ass_pos_cells_place{Ass_2}; 

a_ord = cool(4); 


figure(900)
clf
maximize
ax(1) = subplot(7, 1, 1);
hold on
plot(ms_trk.time/1000, behav.position(:,1));
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylabel('position on track (cm)')
set(gca, 'XTick', []);

ax(2) = subplot(7,1,2:3);
cla
hold on
% assembly 1
for ii = size(Ass_1_ids, 1):-1:1
    s_t = ms_trk.time(Csp(:,Ass_1_ids(ii)) >0)/1000;
    if ~isempty(s_t)
       plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', c_ord(Ass_1,:), 'linewidth', 2)
    end
end
% assembly 2
for ii = size(Ass_2_ids, 1):-1:1
    s_t = ms_trk.time(Csp(:,Ass_2_ids(ii)) >0)/1000;
    if ~isempty(s_t)
       plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5+size(Ass_1_ids, 1), (ones(size(s_t))*ii)+.5+size(Ass_1_ids, 1)]', 'color', c_ord(Ass_2,:), 'linewidth', 2)
    end
end

non_ass_idx = 1:size(gau_sdf, 2); 
rm_idx = (ismember(non_ass_idx, Ass_1_ids)) | (ismember(non_ass_idx, Ass_2_ids)); 
non_ass_idx(Ass_1_ids) = []; 
offset = size(Ass_1_ids, 1) + size(Ass_2_ids, 1);

for ii = size(non_ass_idx, 2):-1:1
    s_t = ms_trk.time(Csp(:,non_ass_idx(ii)) >0)/1000;
    if ~isempty(s_t)
       plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5+offset, (ones(size(s_t))*ii)+.5+offset]', 'color', [.7 .7 .7], 'linewidth', 1)
    end
end
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylim([0 size(Csp, 2)])
set(gca, 'XTick', [], 'YDir', 'normal', 'ytick', [])
ylabel('Cell activity')


ax(3) = subplot(7,1,4:6);
cla
hold on
% assembly 1
for ii = size(Ass_1_ids, 1):-1:1
    
    plot(ms.time/1000, ms_trk_rem.RawTraces(:,Ass_1_ids(ii))*2+ii, 'color', c_ord(Ass_1,:), 'linewidth', 1)
end
% assembly 2
for ii = size(Ass_2_ids, 1):-1:1
    
    plot(ms.time/1000, ms_trk_rem.RawTraces(:,Ass_2_ids(ii))*2+ii+size(Ass_1_ids, 1), 'color', c_ord(Ass_2,:), 'linewidth', 1)
end
set(gca, 'XTick', [], 'YDir', 'normal', 'ytick',[(size(Ass_1_ids, 1)+size(Ass_2_ids, 1))/4, (size(Ass_1_ids, 1)+size(Ass_2_ids, 1)) - (size(Ass_1_ids, 1)+size(Ass_2_ids, 1))/4],...
    'YTickLabel', {['Assembly #' num2str(Ass_1)], ['Assembly #' num2str(Ass_2)]}, 'YTickLabelRotation', 0)
ylabel('Ca2 activity')

ax(4) = subplot(7,1,7);
cla
hold on

    plot(tvec, zscore(time_proj_pos_place(Ass_1,:)), 'color', c_ord(Ass_1,:))

    plot(tvec, zscore(time_proj_pos_place(Ass_2,:)), 'color', c_ord(Ass_2,:))

% view(0,  0)
ylabel({'zscore'; 'react strength'})
xlabel('time (s)')

% view(0,  0)
legend({['Assembly #' num2str(Ass_1)], ['Assembly #' num2str(Ass_2)]})
% set(gca, 'YTick', [0 20], 'YTickLabel', {['Assembly #' num2str(Ass_1)], ['Assembly #' num2str(Ass_2)]}, 'YTickLabelRotation', 0)

% ylabel('assembly ID')
% xlabel('time (s)')
% ylim([min(zscore(time_proj_pos_place(1,:))+10), max(zscore(time_proj_pos_place(1,:))+10*size(time_proj_pos_place,1))])
% xlim([tvec(1) tvec(end)])
%     xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
%     zlim([0 size(time_proj_pos,1)*100])
%     set(gca, 'color', 'k')

% plot the SCE
% ax(4) = subplot(6,1,6);
% cla
% hold on
% plot(ms.time/1000, pop_act, 'color', c_ord(10,:), 'linewidth', 2)
% plot(ms.time/1000, nanmean(all_shuff) - shuff_sd*3,'--',  'color', [.8 .8 .8 .3], 'linewidth',.5)
% plot(ms.time/1000, nanmean(all_shuff) + shuff_sd*3,'--', 'color', [.8 .8 .8 .3], 'linewidth',.5)
% plot(ms.time/1000, nanmean(all_shuff), 'color', [.8 .8 .8 .3], 'linewidth', 1)
% 
% ylabel('pop activity')


linkaxes(ax, 'x')
xlim([ms.time(1)/1000 ms.time(end)/1000])

% save the over all
   exportgraphics(gcf, [data_dir filesep 'Assembly' filesep   sub '_' sess  '_wake.pdf'], 'ContentType', 'vector');

subplot(7,1,2:3)
ylim([0.5 offset+.5])
xlim([0 60])
   exportgraphics(gcf, [data_dir filesep 'Assembly' filesep   sub '_' sess  '_wake_zoom.pdf'], 'ContentType', 'vector');
   
   close(900)
%% Get the wake stem and place plots for the example assemblies. 

for iA = [Ass_1, Ass_2]

    Ass_t_maps = Ass_map{iA}; 

    
l = 1:12; 

m = 1:2:l(end)*2;

n = 2:2:l(end)*2;

figure(iA);clf;
subplot(l(end), 2, m)

hold on
hold on
stem(Ass_pos(:,iA), 'color',[0 0 0])
view(90,90)


p_ord = parula(size(Ass_t_maps,1)); 

for ii = 1:size(Ass_t_maps,1)
    subplot(l(end), 2, m)

    stem(Ass_pos_cells_place{iA}(ii), Ass_pos(Ass_pos_cells_place{iA}(ii),iA),'o','MarkerSize', 15,  'color', p_ord(ii,:), 'Markerfacecolor', p_ord(ii,:))

    
    subplot(l(end), 2, n(ii))
    area(p_bins_int, Ass_t_maps(ii,:), 'facecolor', p_ord(ii,:), 'edgecolor', p_ord(ii,:))
    xlim([p_bins(1) p_bins(end)])
    set(gca, 'xtick', [], 'ytick', [0 1])
    
%     ylim([.5 size(this_ass_map,1)+.5])
end
    subplot(l(end), 2, n(end))
    area(p_bins_int, mean(Ass_t_maps), 'facecolor', c_ord(iA,:), 'edgecolor', c_ord(iA,:))
    xlim([p_bins(1) p_bins(end)])
    ylim([0 1])
    set(gca, 'xtick', [], 'ytick', [0 1])

    
    subplot(l(end), 2, m)

stem(Ass_pos(:,iA), 'color', c_ord(iA,:))
view(90,90)


stem(Ass_pos_cells_place{iA}, Ass_pos(Ass_pos_cells_place{iA},iA), 'color', c_ord(iA,:), 'MarkerFaceColor', c_ord(iA,:))
xlim([0 length(Ass_pos(:,iA))])
set(gca, 'XDir', 'normal')
title(['Assembly #' num2str(iA)])
        

   exportgraphics(gcf, [data_dir filesep 'Assembly' filesep   sub '_' sess  '_Asb_' num2str(iA) '.pdf'], 'ContentType', 'vector');
% close(iA)
end

%% REM reactivation plot for example assemblies and overall 


figure(999)
clf
maximize
ax(1) = subplot(7, 1, 1);
cla
imagesc(rem_time, decode_bin,   decode_prob)
ylabel('decoded position')
set(gca, 'xtick', [], 'YDir', 'normal')
c = colorbar('Location', 'east');
c.AxisLocation = 'out';
pos = c.Position;
c.Position = [pos(1)+.025 pos(2) pos(3)/1.5 pos(4)];
c.Label.String = {'decoded'; 'probability'};
c.Label.FontSize = 14;
caxis([0 .3])
title('Reactivation of wake assemblies in REM')


ax(2) = subplot(7,1,2:3);
cla
hold on
% assembly 1
for ii = size(Ass_1_ids, 1):-1:1
    s_t = rem_time(all_binary_post_REM(:,Ass_1_ids(ii)) >0);
    if ~isempty(s_t)
       plot([s_t; s_t], [(ones(size(s_t))*ii)-.5; (ones(size(s_t))*ii)+.5], 'color', c_ord(Ass_1,:), 'linewidth', 1)
    end
end
% assembly 2
for ii = size(Ass_2_ids, 1):-1:1
    s_t = rem_time(all_binary_post_REM(:,Ass_2_ids(ii)) >0);
    if ~isempty(s_t)
       plot([s_t; s_t], [(ones(size(s_t))*ii)-.5+size(Ass_1_ids, 1); (ones(size(s_t))*ii)+.5+size(Ass_1_ids, 1)], 'color', c_ord(Ass_2,:), 'linewidth', 2)
    end
end

non_ass_idx = 1:size(all_binary_post_REM, 2); 
rm_idx = (ismember(non_ass_idx, Ass_1_ids)) | (ismember(non_ass_idx, Ass_2_ids)); 
non_ass_idx(Ass_1_ids) = []; 
offset = size(Ass_1_ids, 1) + size(Ass_2_ids, 1);

for ii = size(non_ass_idx, 2):-1:1
    s_t = rem_time(all_binary_post_REM(:,non_ass_idx(ii)) >0);
    if ~isempty(s_t)
       plot([s_t; s_t], [(ones(size(s_t))*ii)-.5+offset; (ones(size(s_t))*ii)+.5+offset], 'color', [.7 .7 .7], 'linewidth', 1)
    end
end
xlim([rem_time(1) rem_time(end)])
ylim([0 size(Csp_rem, 2)])
set(gca, 'XTick', [], 'YDir', 'normal', 'ytick', [])
ylabel('Cell activity')


ax(3) = subplot(7,1,4:6);
cla
hold on
% assembly 1
for ii = size(Ass_1_ids, 1):-1:1
    
    plot(rem_time, all_RawTraces_post_REM(:,Ass_1_ids(ii))*2+ii, 'color', c_ord(Ass_1,:), 'linewidth', 1)
end
% assembly 2
for ii = size(Ass_2_ids, 1):-1:1
    
    plot(rem_time, all_RawTraces_post_REM(:,Ass_2_ids(ii))*2+ii+size(Ass_1_ids, 1), 'color', c_ord(Ass_2,:), 'linewidth', 1)
end
set(gca, 'XTick', [], 'YDir', 'normal', 'ytick',[(size(Ass_1_ids, 1)+size(Ass_2_ids, 1))/4, (size(Ass_1_ids, 1)+size(Ass_2_ids, 1)) - (size(Ass_1_ids, 1)+size(Ass_2_ids, 1))/4],...
    'YTickLabel', {['Assembly #' num2str(Ass_1)], ['Assembly #' num2str(Ass_2)]}, 'YTickLabelRotation', 0)
ylabel('Ca2 activity')

ylim([0 (size(Ass_1_ids, 1)+size(Ass_2_ids, 1))+2])
xlim([rem_time(1) rem_time(end)])


ax(4) = subplot(7,1,7);
cla
hold on

    this_proj = (wake_time_proj_rem(Ass_1,:) - mean(time_proj_pos_place(Ass_1,:)))/std(time_proj_pos_place(Ass_1,:)); 
    plot(tvec_rem, this_proj, 'color', c_ord(Ass_1,:))
    
    this_proj = (wake_time_proj_rem(Ass_2,:) - mean(time_proj_pos_place(Ass_2,:)))/std(time_proj_pos_place(Ass_2,:)); 
    plot(tvec_rem, this_proj, 'color', c_ord(Ass_2,:))
    
ylabel({'zscore (to wake)'; 'react strength'})
xlabel('time (s)')

% view(0,  0)
legend({['Assembly #' num2str(Ass_1)], ['Assembly #' num2str(Ass_2)]})




linkaxes(ax, 'x')
xlim([rem_time(1) rem_time(end)])

% save the over all
   exportgraphics(gcf, [data_dir filesep 'Assembly' filesep   sub '_' sess  '_REM.pdf'], 'ContentType', 'vector');

subplot(7,1,2:3)
ylim([0.5 offset+.5])
[~, idx] = max(this_proj); 
xlim([rem_time(idx) - 30 rem_time(idx)+30])
   exportgraphics(gcf, [data_dir filesep 'Assembly' filesep   sub '_' sess  '_REM_zoom.pdf'], 'ContentType', 'vector');
   
   close(999)

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