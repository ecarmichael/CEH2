%% sandbox_PCA/ICA

codebase_dir = 'C:\Users\ecarm\Documents\GitHub\'; 
code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';

RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods'; 

data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3'; 
restoredefaultpath
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));
addpath(code_dir)

cd(data_dir)

%% 

load('ms_trk.mat')

%%
Csp = ms_trk.deconv./ms_trk.denoise; 
Csp = Csp > 0.1; 
ms_trk.Csp = Csp; 



cfg_plot.Ca_type = 'Csp'; 
MS_plot_ca(cfg_plot, ms_trk)

%% bin and convolve
binsize = 0.1; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window

gau_z = zscore(gau_sdf, [], 2); 
%% plot the gaussian smoothed Csp
figure(101)
imagesc(ms_trk.time, 1:size(Csp,2), gau_z')



%% try the assembly code....

Ass_Temp = assembly_patterns(gau_z);

time_proj = assembly_activity(Ass_Temp,gau_sdf'); 

%% plot the output

figure(202)
ax(1) = subplot(5,1,1:4);
hold on
for ii = size(gau_sdf, 2):-1:1
    s_t = ms_trk.time(Csp(:,ii) >0);
    if ~isempty(s_t)
       plot([s_t, s_t]', [(ones(size(s_t))*ii)-.5, (ones(size(s_t))*ii)+.5]', 'color', 'k', 'linewidth', 4)
        
    end
%     plot(ms_trk.time, gau_sdf(:,ii)+ii)
end
xlim([ms_trk.time(1) ms_trk.time(end)])
ylim([0 size(Csp, 2)])
    
ax(2) = subplot(5,1,5);
hold on
c_ord = linspecer(10); 
for ii = 1:10
   plot(ms_trk.time, time_projection(ii,:), 'color', c_ord(ii,:)) 
    
end
    xlim([ms_trk.time(1) ms_trk.time(end)])

linkaxes(ax, 'x')





