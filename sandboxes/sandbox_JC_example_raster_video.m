addpath(genpath('/home/williamslab/Documents/Github/CEH2'));


inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';

% pv1060 LTD5
raw_dir = '/home/williamslab/Desktop/7_19_2019_PV1060_LTD5/H13_M18_S32_LTD5';
decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/LTD5';
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/LTD5';


% raw_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1060/11_23_2019_PV1060_HATD5/H13_M25_S1_HATD5';
% decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/HATD5';
% % ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
% ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/HATD5';
% 


core_colors = linspecer(4); 
%% get the decoder files

cd(decode_dir)
load('decoding.mat')
if ~isfield(decoding, 'WAKE_decoded_position')
    decoding.WAKE_decoded_position = [];
    [max_prob, pos_idx] = max(decoding.WAKE_decoded_probabilities,[], 1);
    for ii = length(pos_idx):-1:1
        if isnan(max_prob(ii))
            decoding.WAKE_decoded_position(ii) = NaN;
        else
            decoding.WAKE_decoded_position(ii) = decoding.bin_centers_vector(pos_idx(ii));
        end
    end
end

%%  Get the raw videos
% clearvars -except raw_dir ms_dir
cd(raw_dir)

vid_num = 2;

b_obj = VideoReader(['behavCam' num2str(vid_num) '.avi']);
for iF = b_obj.NumFrames:-1:1
    b_f{iF} = read(b_obj, iF);
    b_t(iF) = b_obj.CurrentTime;
end

ca_obj = VideoReader(['msCam' num2str(vid_num) '.avi']);
for iF = ca_obj.NumFrames:-1:1
    ca_f{iF} = read(ca_obj, iF);
    ca_t(iF) = ca_obj.CurrentTime;
end

%%
% cd(ms_dir)
warning off; load('ms.mat'); warning on; 
ms = msExtractBinary_detrendTraces(ms);
% frame_offset = ms.timestamps(1);
Fs = mode(diff(ms.time));

vid_idx = [1 find(diff(ms.frameNum) ~=1)+1]; % get the indicies where the video starts.
this_idx = vid_idx(vid_num);
if vid_num == length(vid_idx)
    next_idx = length(ms.RawTraces);
else
    next_idx = vid_idx(vid_num+1);
end

%% limited cell version
% good_cells =20;
% % find x nice cells.
% if exist('SA.mat', 'file')
%     load('SA.mat')
%     cells_to_use = SA.WholePlaceCell;
%     c_ord = linspecer(good_cells);
% elseif exist('with modfied eva/NonPlaceCell.mat', 'file')
%     load('with modfied eva/NonPlaceCell.mat')
%     %     for iS = length(spatial_analysis.bin):-1:1
%     %        cells_to_use(iS) = spatial_analysis.bin{iS, 3}.IsPlaceCell;
%     %     end
%     cells_to_use = ones(size(ms.RawTraces,2),1);
%     cells_to_use(NonPlaceCell.CellID) = 0;
%     cells_to_use = find(logical(cells_to_use));
%     c_ord = linspecer(good_cells+6);
%     c_ord(((good_cells/2)-1):(good_cells/2)+4,:) = []; % remove yellows. They look terrible.
% else
%     cells_to_use = 1:size(ms.RawTraces,2);
%     c_ord = linspecer(good_cells+6);
%     c_ord(((good_cells/2)-1):(good_cells/2)+4,:) = []; % remove yellows. They look terrible.
% end
% 
% d_max = []; 
% for iC = length(cells_to_use):-1:1
%        this_data = zscore(ms.RawTraces(this_idx:next_idx,cells_to_use(iC)));
%     pks = findpeaks(zscore(ms.RawTraces(this_idx:next_idx,cells_to_use(iC))),33,'MinPeakDistance',2, 'MinPeakHeight', 5);
%         d_max(iC) = max(this_data);
%     d_pks(iC) = length(pks);
% end
% [best_pks, best_cells] = sort(d_pks, 'descend');
% [best_most_pks, best_most_cells] = sort(d_pks(best_cells(1:40)), 'descend');
% cells = best_cells(1:good_cells);
% % cells = 1:good_cells;
% 
% sub_mat = reshape(1:((6+length(cells))*4),4,6+length(cells))'; % matrix to pull subplot values from
%%
if ishandle(108)
    close(108)
end
figure(108)
subplot(6,4,[1:4])
    set(gca, 'color', 'k'); axis off; box on;
hold on
imagesc(b_f{1})
if contains(raw_dir, 'HAT')
    rectangle('position', [0 0 50 1], 'facecolor', [core_colors(2,:) .4], 'EdgeColor', [core_colors(3,:), 0]);
    rectangle('position', [50 0 50 1], 'facecolor', [core_colors(3,:) .4], 'EdgeColor', [core_colors(3,:), 0])
end
title(sprintf('Time = %.2fs', ms.time(1)/1000));
set(gca, 'xtick', [], 'ytick', [])

% plot the decoder
% ax(100) = subplot(6,4,5:6);
% imagesc(ms.time/1000, 1:size(decoding.WAKE_decoded_probabilities,1), decoding.WAKE_decoded_probabilities)
% % caxis([0 max(decoding.REM_decoded_probabilities(3:end-3,:), [], 'all')])

set(gcf, 'color', 'k')

ax(200)= subplot(6,4,[ 5 6 9 10 13 14 17 18 21 22]);
MS_Ca_Raster(ms.Binary(this_idx:next_idx,:)', ms.time(this_idx:next_idx)/1000)
linkaxes(ax, 'x');
set(gca, 'color', 'k', 'XColor',[.6 .6 .6], 'YColor',[.6 .6 .6] ); %set background color. 
ylabel('Cell ID', 'Color', [.6 .6 .6])
xlabel('time (ms)', 'Color', [.6 .6 .6])

subplot(6,4,[7 8 11 12 15 16 19 20 23 24])
imagesc(ca_f{1})
colormap('gray')
set(gca, 'xtick', [], 'ytick', [])
hold on
% for iC = length(cells):-1:1
%     [x(iC),y(iC)]=find(ms.SFPs(:,:,cells(iC)) == max(ms.SFPs(:,:,cells(iC)),[],'all'));
% end
% 
% scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5)

%%
set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')
max_pos = max(decoding.WAKE_decoded_position); 
max_deco_prob = max(decoding.WAKE_decoded_probabilities,[], 'all');

clear F
for iF = Fs:1:size(b_t,2)-Fs
    
    subplot(6,4,1:4)
    hh = imagesc(b_f{iF});
    title(sprintf('Time = %.3fs', ms.time(this_idx+iF)/1000));
    set(gca, 'xtick', [], 'ytick', [])
    hold on
    if isnan(decoding.WAKE_decoded_position(iF))
         h = [];
    else
        h = scatter((decoding.WAKE_decoded_position(this_idx+iF)/max_pos)*(size(b_f{iF},2)), size(b_f{iF},1)/2, 800, 'LineWidth', 5);
        h.MarkerEdgeColor = core_colors(2,:);
%         h.MarkerEdgeAlpha  =max(decoding.WAKE_decoded_probabilities(:,iF))/max_deco_prob;
    end
    %     title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(iF), ms.frameNum(this_idx+iF)*Fs));
    subplot(6,4,[5 6 9 10 13 14 17 18 21 22])
    xlim([ms.time(this_idx+iF)/1000 - 2.5 ms.time(this_idx+iF)/1000 + 2.5])
        hl = vline(median([ms.time(this_idx+iF)/1000 - 2.5 ms.time(this_idx+iF)/1000 + 2.5]));
    
    subplot(6,4,[7 8 11 12 15 16 19 20 23 24])
    imagesc(ca_f{iF})
%     hold on
%     scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5);
    set(gca, 'xtick', [], 'ytick', [])
    
    drawnow;
    F(iF) = getframe(gcf) ;
    delete(hh);
    delete(h); 
    delete(hl)
    
end


% set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

%%
parts = strsplit(cd, filesep);
mkdir(['/home/williamslab/Dropbox (Williams Lab)/Decoding_data/Videos/' parts{end-1}])
writerObj = VideoWriter(['/home/williamslab/Dropbox (Williams Lab)/Decoding_data/Videos/' parts{end-1} filesep 'Example_Ca_WAKE_decode_raster.avi']);
rate = 1;
writerObj.FrameRate = Fs/rate;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = Fs:rate:size(b_t,2)-Fs
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


