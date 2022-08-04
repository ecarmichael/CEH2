%% sandbox_plot_CA_trace_JC
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared')); % where the codebase repo can be found
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2')); % where the multisite repo can be found

data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATD1';
cd(data_dir); 

load('ms_resize.mat')

load('ms_trk.mat')

% load('all_RawTraces_pre_REM.mat')

% load('all_RawTraces_post_REM.mat')
ms_seg_resize.hypno_label
%% split into pre, track, and post
close all
cfg_plot = [];
cfg_plot.plot_type = '2d'; 
% cfg_plot.rescale = 'max';
cfg_plot.Ca_chan = [14 23 34 59 69 76 77 88 100 111 122 125 129 131 146 147 157 162 187 200 229 243 245 262 273];
% cfg_plot.Ca_chan = cfg_plot.Ca_chan+25
% cfg_plot.Ca_chan = 250:300;
cfg_plot.offset = .7;
% cfg_plot.colors = repmat([0 0 0], length(cfg_plot.Ca_chan),1)
cfg_plot.colors = winter(length(cfg_plot.Ca_chan)+5);

ax(1) = subplot(1,3,1);
temp_ms = [];
temp_ms.time = ms_seg_resize.time{7};
temp_ms.RawTraces = ms_seg_resize.RawTraces{7}; 
MS_plot_ca(cfg_plot, temp_ms)
xlabel('time (s)');
ylabel('cell id')
yt = get(gca, 'ytick'); yl = get(gca, 'yticklabel');
set(gca, 'ytick', [yt(1) yt(end)])
set(gca, 'yticklabel', {yl{1}; yl{end}})
xlim([0 60])

ax(2) = subplot(1,3,2);
MS_plot_ca(cfg_plot, ms_trk)
yt = get(gca, 'ytick'); yl = get(gca, 'yticklabel');
% set(gca, 'ytick', [yt(1) yt(end)])
% set(gca, 'yticklabel', {yl{1}; yl{end}})
xlim([0 60])


ax(3) = subplot(1,3,3);
temp_ms = [];
temp_ms.time = ms_seg_resize.time{15};
temp_ms.RawTraces = ms_seg_resize.RawTraces{15}; 
MS_plot_ca(cfg_plot, temp_ms)
yt = get(gca, 'ytick'); yl = get(gca, 'yticklabel');
% set(gca, 'ytick', [yt(1) yt(end)])
% set(gca, 'yticklabel', {yl{1}; yl{end}})
xlim([0 60])

linkaxes(ax, 'y');
% ylim([0 ((cfg_plot.Ca_chan(end) - cfg_plot.Ca_chan(1))*cfg_plot.offset)*1.1])

%%
saveas(gcf, ['C:\Users\ecarm\Desktop' filesep 'pv160_HATD1_trace'], 'png');
saveas(gcf, ['C:\Users\ecarm\Desktop' filesep 'pv160_HATD1_trace'], 'fig');
saveas(gcf, ['C:\Users\ecarm\Desktop' filesep 'pv160_HATD1_trace'], 'svg');