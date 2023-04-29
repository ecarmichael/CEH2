%% sandbox_PCA/ICA

codebase_dir = 'C:\Users\ecarm\Documents\GitHub\'; 
code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';

RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods'; 

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3'; 
restoredefaultpath
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(code_dir)

cd(data_dir)

%% 

load('ms_trk.mat')
load('behav_DLC.mat')

% remove inactive cells
keep_idx = sum(ms_trk.deconv, 1) >0; 

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv'; 
ms_trk_rem = MS_Remove_trace(cfg_rem, ms_trk); 


%% follow grosmark et al. method of deconv preprocessing 
Csp = ms_trk_rem.deconv./ms_trk_rem.denoise; 
Csp = Csp > 0.01; 
ms_trk_rem.Csp = Csp; 

cfg_plot.Ca_type = 'Csp'; 
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

%% try the assembly code....

Ass_Temp = assembly_patterns(gau_z);

time_proj = assembly_activity(Ass_Temp,gau_sdf'); 


%% color code the assemblies 
Ass_sort = []; figure(303); hold on
for ii = size(Ass_Temp,2):-1:1
    
    [Ass_sort(:,ii), this_idx] = sort(Ass_Temp(:,ii), 'descend'); 
    thresh = prctile(Ass_Temp(:,ii), 95); 
    c_idx = Ass_Temp(this_idx,ii) > thresh; 
    
%     plot(Ass_Temp(c_idx, ii), )
    
    
end

%% stem plot for first few ensembles 
 figure(303); hold on
for ii = 2
        stem(Ass_Temp(:,ii))
end

%% plot the output

figure(202)
ax(1) = subplot(6, 1, 1)
hold on
plot(ms_trk.time/1000, behav.position(:,1)); 
plot(ms_trk.time/1000, behav.position(:,2)); 
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ax(2) = subplot(6,1,2:5);
hold on
for ii = size(gau_sdf, 2):-1:1
    s_t = ms_trk.time(Csp(:,ii) >0)/1000;
    if ~isempty(s_t)
       plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', 'k', 'linewidth', 4)
        
    end
%     plot(ms_trk.time, gau_sdf(:,ii)+ii)
end
xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])
ylim([0 size(Csp, 2)])
    
ax(3) = subplot(6,1,6);
hold on
c_ord = linspecer(10); 
for ii = 1:10
   plot3(ms_trk.time/1000, zscore(time_proj(ii,:)),ones(1, length(time_proj))*100* ii, 'color', c_ord(ii,:)) 
    
end
    xlim([ms_trk.time(1)/1000 ms_trk.time(end)/1000])

linkaxes(ax, 'x')





