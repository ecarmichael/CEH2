%% sandbox_PCA/ICA

codebase_dir = 'C:\Users\ecarm\Documents\GitHub\'; 
code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';

RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods'; 

% data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3'; 
data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter'
restoredefaultpath
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(code_dir)

cd(data_dir)

%% 

% load('ms_trk.mat')
% load('behav_DLC.mat')

load('pv1069_LTD5_data.mat')

behav = MS_align_data(behav,ms)

ms_trk = ms; 
keep_idx = zeros(1,size(ms_trk.RawTraces,2));
keep_idx(1:floor(size(ms_trk.RawTraces,2)*.66)) = 1; 


cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'RawTraces'; 
ms_trk_rem = MS_Remove_trace(cfg_rem, ms_trk); 


ms_trk_rem = MS_append_deconv(ms_trk_rem); 

% remove inactive cells

keep_idx = sum(ms_trk_rem.deconv, 1) >0; 

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv'; 
ms_trk_rem = MS_Remove_trace(cfg_rem, ms_trk_rem); 

% deconv the all_b



%% follow grosmark et al. method of deconv preprocessing 
Csp = ms_trk_rem.deconv./ms_trk_rem.denoise; 
Csp = Csp > 0.01; 
ms_trk_rem.Csp = Csp; 

cfg_plot.Ca_type = 'RawTraces'; 
cfg_plot.plot_type = '2d'; 
MS_plot_ca(cfg_plot, ms_trk_rem)

%% bin and convolve
binsize = 0.1; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window

gau_z = zscore(gau_sdf, [], 2); 
%% plot the gaussian smoothed Csp
figure(101)
ax(1)= subplot(5,1,1);
hold on
plot(ms_trk.time, behav.position(:,1)); 
plot(ms_trk.time, behav.position(:,2)); 
xlim([ms_trk.time(1) ms_trk.time(end)])

ax(2) = subplot(5,1,2:5);
imagesc(ms_trk.time, 1:size(Csp,2), gau_z')

linkaxes(ax, 'x')


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

this_cell = ms.time(find(Csp(:,ii)))/1000;

spk_count = histc(this_cell,tbin_edges); % get spike counts for each bin
spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.

data_h(:,ii) = spk_count; 
end



tvec = tbin_centers; 


%% convert data if needed. 
gau_z_t = gau_z;
gau_sdf_t = gau_sdf;

gau_z = data_int; 
gau_sdf = data_int; 

%% try the assembly code....

Ass_Temp = assembly_patterns(data_h');

time_proj = assembly_activity(Ass_Temp,data_h'); 


%% color code the assemblies 
Ass_sort = []; 
%figure(303);clf;  hold on
for ii = size(Ass_Temp,2):-1:1
    
    [Ass_sort(:,ii), this_idx] = sort(Ass_Temp(:,ii), 'descend'); 
    thresh = prctile(Ass_Temp(:,ii), 95); 
    c_idx = Ass_Temp(this_idx,ii) > thresh; 
    
%     plot(Ass_Temp(c_idx, ii))
%             stem(Ass_Temp(c_idx, ii))

end

%% keep assemblies with strong activations
keep_idx = zeros(1,size(Ass_Temp,2));  
for ii = size(Ass_Temp,2):-1:1

   if max(Ass_Temp(:, ii)) > 0.2
       keep_idx(ii) = 1; 
   end
    
end

Ass_pos = Ass_Temp; 
Ass_pos(:,~keep_idx) = []; 

time_proj_pos = time_proj;
time_proj_pos(~keep_idx,:) = [];
%% stem plot for first few ensembles 
 figure(303); hold on
 c_ord = parula(8);
for ii = 1:8
    subplot(4,2,ii)
        stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
        view(90,90)
end

%% plot the output

figure(202)
clf
% maximize
ax(1) = subplot(6, 1, 1)
hold on
plot(ms_trk.time/1000, behav.position(:,1)); 
% plot(ms_trk.time/1000, behav.position(:,2)); 
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylabel('position on track (cm)')
set(gca, 'XTick', []);

ax(2) = subplot(6,1,2:4);
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

ax(3) = subplot(6,1,5:6);
cla
hold on
c_ord = cool(size(time_proj_pos,1)+4); 
for ii = 1:size(time_proj_pos,1)
%    nan_idx = time_proj_pos(ii,:) < prctile(time_proj_pos(ii,:), 10);
%    plot3(ms_trk.time(~nan_idx)/1000, zscore(time_proj_pos(ii,~nan_idx)),ones(1, sum(~nan_idx))*100* ii, '.', 'color', c_ord(ii,:)) 
%        plot3(tvec, zscore(time_proj_pos(ii,:)),ones(1, length(time_proj_pos(ii,:)))*100* ii, 'color', c_ord(ii,:)) 
       plot(tvec, zscore(time_proj_pos(ii,:))+ii*10, 'color', c_ord(ii+4,:)) 

       tick_val(ii) = mode(zscore(time_proj_pos(ii,:))+ii*10); 
       tick_label{ii} = num2str(ii); 
end
% view(0,  15)
set(gca, 'YTick', tick_val, 'YTickLabel', tick_label)

ylabel('assembly ID')
xlabel('time (s)')
ylim([min(zscore(time_proj_pos(1,:))+10), max(zscore(time_proj_pos(1,:))+10*size(time_proj_pos,1))])
xlim([tvec(1) tvec(end)])
%     xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
%     zlim([0 size(time_proj_pos,1)*100])
%     set(gca, 'color', 'k')

linkaxes(ax, 'x')




%% get the assembly triggered position average. 
win = floor(2.5 * mode(diff(behav.time))); 
figure(5)
clf
n = ceil(size(time_proj_pos,1)/2); 
m = 2; 

for ii = 1: size(time_proj_pos,1)
    subplot(m,n,ii)
    hold on
    [~, p_idx] = findpeaks(zscore(time_proj_pos(ii,:)),'MinPeakHeight', 1.96   ,'MinPeakDistance', 2*floor(mode(diff(ms.time)))); 
    
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

%% final figure with coloured assemblies for # 3 and 9. 
Ass_idx = [3 9]; 

figure(202)
clf
% maximize
ax(1) = subplot(6, 1, 1)
hold on
plot(ms_trk.time/1000, behav.position(:,1)); 
% plot(ms_trk.time/1000, behav.position(:,2)); 
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylabel('position on track (cm)')
set(gca, 'XTick', []);

ax(2) = subplot(6,1,2:4);
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

ax(3) = subplot(6,1,5:6);
cla
hold on
c_ord = cool(size(time_proj_pos,1)+4); 
for ii = 1:size(time_proj_pos,1)
%    nan_idx = time_proj_pos(ii,:) < prctile(time_proj_pos(ii,:), 10);
%    plot3(ms_trk.time(~nan_idx)/1000, zscore(time_proj_pos(ii,~nan_idx)),ones(1, sum(~nan_idx))*100* ii, '.', 'color', c_ord(ii,:)) 
%        plot3(tvec, zscore(time_proj_pos(ii,:)),ones(1, length(time_proj_pos(ii,:)))*100* ii, 'color', c_ord(ii,:)) 
       plot(tvec, zscore(time_proj_pos(ii,:))+ii*10, 'color', c_ord(ii+4,:)) 

       tick_val(ii) = mode(zscore(time_proj_pos(ii,:))+ii*10); 
       tick_label{ii} = num2str(ii); 
end
% view(0,  15)
set(gca, 'YTick', tick_val, 'YTickLabel', tick_label)

ylabel('assembly ID')
xlabel('time (s)')
ylim([min(zscore(time_proj_pos(1,:))+10), max(zscore(time_proj_pos(1,:))+10*size(time_proj_pos,1))])
xlim([tvec(1) tvec(end)])
%     xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
%     zlim([0 size(time_proj_pos,1)*100])
%     set(gca, 'color', 'k')

linkaxes(ax, 'x')

