
clear all
cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Decoding/pv1060/HATDSwitch/')
load('1000shuffling.mat', 'Final_start_frame')
load('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/4.PlaceCell/pv1060/HATDSwitch/spatial_analysis.mat')
load('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/JisooProject2020/2020_Results_aftercutting/Across_episodes/Inter/PV1060/11_26_2019_PV1060_HATSwitch/all_binary_post_REM.mat')

parts = strsplit(cd, filesep);
sub = parts{end-1};
sess = parts{end}; 

if exist('spatial_analysis', 'var')
    place_idx = zeros(length(spatial_analysis.bin),1); % allocate the index array
    centroids = nan(size(place_idx));
    
    for iC = length(spatial_analysis.raw):-1:1
        if spatial_analysis.bin{iC,3}.IsPlaceCell
            place_idx(iC) = 1;
            centroids(iC) = spatial_analysis.bin{iC,3}.PlaceFieldCentroid{1}(1);
        end
    end
    centroids = (centroids*4)+10; 
   [~, cent_sort] = sort(centroids); 
   place_idx = place_idx(cent_sort); 

else
    cent_sort = 1:size(all_binary_post_REM,2); % just use default sort. 
    place_idx = ones(size(cent_sort)); 
end

win_s = [1,1]; 

%%
for ii = 1:length(Final_start_frame)

% make a raster
this_ca = all_binary_post_REM(Final_start_frame(ii) - win_s(1)*30:Final_start_frame(ii) + win_s(2)*30,:)';
this_tvec = Final_start_frame(ii) - win_s(1)*30:Final_start_frame(ii) + win_s(2)*30; 
this_tvec = ((this_tvec - this_tvec(1))/30) - win_s(1);

figure(ii);
c_mat = [linspecer(sum(place_idx));  repmat([1 1 1], sum(~place_idx),1)]; % make colors depending on the
MS_Ca_Raster(this_ca(cent_sort,:),this_tvec, 14, c_mat);%repmat([1,1,1], size(this_ca,1),1)
xlabel('time (s)')
ylabel('Cell ID')
xline(0, '--w', 'start', 'linewidth', 2);
x_lim = xlim;
% xline(x_lim(2) - win_s, '--w', 'end', 'linewidth', 2);
set(gca, 'color', 'k'); %set background color.
colormap([linspecer(sum(place_idx));  repmat([1 1 1], 1,1)]);
cx = colorbar;
caxis([10 90])
if exist('centroids', 'var')
%     cx.Ticks = 10:20:90;
    cx.Label.String = 'place cell centroid';
end

cfg_plot = [];
cfg_plot.ft_size = 12; 
SetFigure(cfg_plot, gcf)
 set(0, 'DefaultFigureRenderer', 'painters');
saveas(gcf, ['/home/williamslab/Desktop/REM_Replay_Rasters' filesep sub '_' sess 'pv1060_HATDS_' num2str(ii)], 'png');

saveas(gcf, ['/home/williamslab/Desktop/REM_Replay_Rasters' filesep sub '_' sess 'pv1060_HATDS_' num2str(ii)], 'svg');
close all
end