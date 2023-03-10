function [out]  = MS_get_SCE_JC(ms, behav, save_dir); 

%% sandbox for SCEs based on malvache

% In order to reduce noise, a third-order Savitzky-Golay filter with a frame size of 500ms was first
% applied on the fluorescence calcium signal of each cell. The threshold for detecting calcium transients
% was adapted for each time point and each cell as follows: it was the sum of the median value with
% three times the interquartile range calculated within a -2/+2s sliding window. To avoid detecting
% twice the same calcium transient, the minimal delay between events was set to one second. Activity
% occurring during run epochs was not included in this analysis. SCE corresponded to the synchronous
% calcium events that involved more cells than expected by chance within a 200ms time window (i.e.
% >3 std after temporal reshuffling of cell activity) and with a minimum 5 cells participation threshold.

% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\11_19_2019_PV1060_HATD1'
% cd(data_dir)
% load('ms_trk.mat')

% behav_dir = 'K:\Jisoo_Project\RawData\pv1060\11_26_2019_PV1060_HATSwitch\H13_M5_S15_HATSwtich';
% cd(behav_dir)
% load('behav.mat'); 
% load('ms.mat')

if nargin < 3
    save_dir = cd;
end

if ~isfield(ms, 'Binary')
ms = msExtractBinary_detrendTraces(ms); 
end


% interp behaviour
behav_i = behav; 
behav_i.time = behav.time - behav.time(1); 
behav_i.time = behav_i.time+ms.time(1); 
behav_i.position = []; behav_i.speed = []; 
[~, u_idx] = unique(behav_i.time); 
behav_i.position(:,1) = interp1(behav_i.time(u_idx), behav.position(u_idx,1), ms.time); 
behav_i.speed = interp1(behav_i.time(u_idx), behav.speed(u_idx), ms.time); 

behav_i.time = ms.time; 


behav_i.speed(behav_i.speed > 70) = NaN; 
behav_i.speed = fillmissing(behav_i.speed, 'spline'); 
behav_i.speed = smooth(behav_i.speed, mode(diff(behav_i.time))); 
behav_i.speed(behav_i.speed < 0) = 0; 


move_idx = behav_i.speed > 10; 

df = mode(diff(ms.time)); 
frame_n = floor(df/2);
if mod(frame_n, 2)==0; frame_n = frame_n+1; end
frame_n_200 = floor(df/5); 

% %% try deconv instead
% 
%  for iChan = size(ms_trk.RawTraces,2):-1:1
%             tic;
%             [denoise,deconv] = deconvolveCa(ms_trk.detrendRaw(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
%             toc;
%             all_denoise(:,iChan) = denoise;    all_deconv(:,iChan) = deconv;
%  end
%     ms_trk.denoise = all_denoise;
%     ms_trk.deconv = all_deconv;
%     ms_trk.decon_bin = all_deconv;
%     ms_trk.decon_bin(all_deconv >0) = 1; 
% %% smooth with S-G filter 
% smooth_data = [];
% 
% for ii = size(data, 2):-1:1
%     smooth_data(:,ii) = sgolayfilt(data(:,ii),3,frame_n);
% 
% %     inter quantile threshold (unclear what median value 
%     
% end

%% select data
data_in =ms.Binary;
data_og =ms.Binary;

data_in(move_idx,:) = NaN;

%% SCE check  >3sd shuffled summed activity 

% d_f_f = []; 
% for ii = size(data_in, 2):-1:1
%     
%     d_f_f(:,ii) = zscore(data_in(:,ii)); 
% end

shuff = 100;
all_shuff = [];
tic
for iS = shuff:-1:1
    this_data = []; 
    for ii = size(data_in, 2):-1:1
        this_data(ii,:) = circshift(data_in(:,ii),floor(MS_randn_range(1,1,1,length(data_in(:,ii)))));
    end % end cells

    all_shuff(iS, :) = movmean(nansum(this_data,1), frame_n_200); 
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

% figure(999)
% ax(1) =subplot(8,2,1:8);
% MS_Ca_Raster(data_og', ms.time/1000);
% % xlabel('frame number')
% ylabel('cell id')
% set(gca, 'xtick', [])
% 
% ax(2) =subplot(8,2,9:10);
% hold on
% plot(ms.time/1000, pop_act);
% plot(ms.time(SCE_idx)/1000, pop_act(SCE_idx), 'x')
% % xlim([ms_trk.time(1) ms_trk.time(end)])
% hline(thresh)
% ylabel('pop activity')
% set(gca, 'xtick', [])
% 
% % plot the speed w/ threshold
% ax(3) =subplot(8,2,11:12);
% plot(ms.time/1000, behav_i.speed);
% hline(7.5);
% ylabel('speed (cm/s)')
% set(gca, 'xtick', [])
% 
% ax(4) =subplot(8,2,13:14);
% hold on
% plot(ms.time/1000, behav_i.position(:,1));
% plot(ms.time(SCE_idx)/1000, behav_i.position(SCE_idx,1), 'xr')
% ylabel({'position on'; 'track (cm)'})
% xlabel('time (s)')
% set(gca, 'xtick', [])
% 
% linkaxes(ax, 'x')
% xlim([ms.time(1)/1000 ms.time(end)/1000])
% 
% subplot(8,2,15);
% histogram(behav_i.position(~move_idx,1),40); 
% xlim([min(behav_i.position(:,1)) max(behav_i.position(:,1))])
% ylabel('# SCEs');
% xlabel('position on track');
% 
% subplot(8,2,16);
[occ_vec, bins] = hist(behav_i.position(~move_idx,1),40); 
occ_nan = find(isnan(occ_vec)); % get unoccupied bins. Shouldn't be any.
occ_vec(occ_nan) = NaN; 
occ_vec = occ_vec.*(1/mode(diff(behav.time))); % conver to time 

SCE_vec = hist(behav_i.position(SCE_idx,1),40);
SCE_vec(occ_nan) = NaN; 

rate_vec = SCE_vec ./ occ_vec; 
% bar(bins,  rate_vec)
% % ylim([0 5]);
% ylabel('SCEs / s')
% xlim([min(behav_i.position(:,1)) max(behav_i.position(:,1))])
% xlabel('position on track');
% 
% 
% cfg_fig = [];
% cfg_fig.ft_size = 12; 
% 
% SetFigure(cfg_fig, gcf)
% 
% 
parts = strsplit(cd, filesep); 

% mkdir(save_dir)
% saveas(gcf, [save_dir filesep 'SCE_' parts{end-1} '_' parts{end}], 'png')

out.rate_vec = rate_vec;
out.SEC_vec = SCE_vec; 
out.bins = bins; 
out.occ_vec = occ_vec; 
out.SCE_idx = SCE_idx; 
out.sub = parts{end-1};
out.sess = parts{end}; 

% %% 2d version
% 
% ax(1) =subplot(9,2,1:8);
% MS_Ca_Raster(data_og', ms.time/1000);
% % xlabel('frame number')
% ylabel('cell id')
% set(gca, 'xtick', [])
% 
% ax(2) =subplot(9,2,9:10);
% hold on
% plot(ms.time/1000, pop_act);
% plot(ms.time(SCE_idx)/1000, pop_act(SCE_idx), 'x')
% % xlim([ms_trk.time(1) ms_trk.time(end)])
% hline(thresh)
% ylabel('pop activity')
% set(gca, 'xtick', [])
% 
% % plot the speed w/ threshold
% ax(3) =subplot(9,2,11:12);
% plot(ms.time/1000, behav_i.speed);
% hline(7.5);
% ylabel('speed (cm/s)')
% set(gca, 'xtick', [])
% 
% ax(4) =subplot(9,2,13:14);
% hold on
% plot(ms.time/1000, behav_i.position(:,:));
% % vline(ms_trk.time(SCE_idx)/1000, 'r', '-')
% plot(ms.time(SCE_idx)/1000, behav_i.position(SCE_idx,1), 'xr')
% ylabel({'position on'; 'track (cm)'})
% xlabel('time (s)')
% set(gca, 'xtick', [])
% 
% linkaxes(ax, 'x')
% xlim([ms.time(1)/1000 ms.time(end)/1000])
% 
% subplot(9,2,[15 17]);
% hold on
% plot(behav_i.position(:,1), behav_i.position(:,2), '.', 'color', [.6 .6 .6])
% plot(behav_i.position(SCE_idx,1), behav_i.position(SCE_idx,2), '.r')
% 
% xlim([min(behav_i.position(:,1)) max(behav_i.position(:,1))])
% ylabel('# SCEs');
% xlabel('position on track');
% 
% subplot(9,2,[16 18]);
% x_edges = min(behav_i.position(:,1)):3: max(behav_i.position(:,1));
% y_edges = min(behav_i.position(:,2)):3: max(behav_i.position(:,2));
% 
% [occ_vec, bins] = histcn(behav_i.position(~move_idx,:),x_edges, y_edges); 
% occ_nan = find(isnan(occ_vec)); % get unoccupied bins. Shouldn't be any.
% occ_vec(occ_nan) = NaN; 
% occ_vec = occ_vec.*(1/mode(diff(behav.time))); % conver to time 
% 
% SCE_vec = histcn(behav_i.position(SCE_idx,:),x_edges, y_edges);
% SCE_vec(occ_nan) = NaN; 
% 
% rate_vec = SCE_vec ./ occ_vec; 
% imagesc(flipud(rate_vec')); 
% % set(gca, 'YDir', 'reverse')
% % ylim([0 5]);
% ylabel('SCEs / s')
% % xlim([min(behav_i.position(:,1)) max(behav_i.position(:,1))])
% xlabel('position on track');
% 
% 
% cfg_fig = [];
% cfg_fig.ft_size = 12; 
% 
% SetFigure(cfg_fig, gcf)