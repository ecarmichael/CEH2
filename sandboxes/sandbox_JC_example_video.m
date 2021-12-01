addpath(genpath('/home/williamslab/Documents/Github/CEH2'));


inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));


raw_dir = '/home/williamslab/Desktop/7_19_2019_PV1060_LTD5/H13_M18_S32_LTD5';
decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/LTD5';
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/LTD5';

core_colors = linspecer(3); 
%%
% clearvars -except raw_dir ms_dir
cd(raw_dir)

vid_num = 3;

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
load('ms.mat')
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

%% deconv debugger.
% % counter_init(size(ms.RawTraces,2),size(ms.RawTraces,2))
% for iChan = 1
%     %     counter(iChan, size(ms.RawTraces,2))
%     tic;
%     [denoise,deconv] = deconvolveCa(ms.detrendRaw(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
%     toc;
%     ms.denoise(:,iChan) = denoise;    ms.deconv(:,iChan) = deconv;
% end
% 
% % debugging
% % iChan = 900;
% if ishandle(iChan)
%     close(iChan)
% end
% 
% figure(iChan+100)
% hold on
% plot(zscore(ms.detrendRaw(:,iChan))./max(zscore(ms.detrendRaw(:,iChan))), 'k');
% plot(ms.denoise(:,iChan)-.2, 'r');
% plot(((ms.deconv(:,iChan)./max(ms.deconv(:,iChan)))*.1) -.2, 'b');
% plot((ms.Binary(:,iChan)*.1)-.2, 'g');
% % MS_plot_ca_trace(ms.FiltTraces(1:33*60,1:50)')
% legend('Detrend', 'OASIS: Denoised', 'OASIS: deconv', 'Binary', 'orientation', 'horizontal', 'location', 'north')

%% plot stuff
% close all
% figure(101)
% subplot(3,4,1:4)
% imagesc(b_f{1})
% title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(1), ca_t(1)));
% set(gca, 'xtick', [], 'ytick', [])
% 
% subplot(3,4,[5 6 9 10])
% % MS_Ca_Raster(ms.Binary(1:ms.timestamps(2),:)'); %, ms.time(1:ms.timestamps(2))'/1000
% set(gca, 'xtick', [])
% ylabel('cell number')
% 
% 
% subplot(3,4,[ 7 8 11 12])
% imagesc(ca_f{1})
% colormap('gray')
% set(gca, 'xtick', [], 'ytick', [])
% 
% set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')
% 
% clear F
% for iF = Fs:4:size(b_t,2)-Fs
%     
%     subplot(3,4,1:4)
%     imagesc(b_f{iF})
%     title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(iF), ca_t(iF)));
%     
%     subplot(3,4,[5 6 9 10])
%     xlim([ms.timestamps(1)-Fs + iF  ms.timestamps(1)+Fs + iF])
%     h = vline(median([ms.timestamps(1)-Fs + iF  ms.timestamps(1)+Fs + iF]));
%     
%     subplot(3,4,[ 7 8 11 12])
%     imagesc(ca_f{iF})
%     
%     drawnow;
%     F(iF) = getframe(gcf) ;
%     delete(h);
%     
% end

%% limited cell version
good_cells =20;
% find x nice cells.
if exist('SA.mat', 'file')
    load('SA.mat')
    cells_to_use = SA.WholePlaceCell; 
    c_ord = linspecer(good_cells);
else
    cells_to_use = 1:size(ms.RawTraces,2); 
    c_ord = linspecer(good_cells+6);
c_ord(((good_cells/2)-1):(good_cells/2)+4,:) = []; % remove yellows. They look terrible. 
end

d_max = []; 
for iC = length(cells_to_use):-1:1
       this_data = zscore(ms.RawTraces(this_idx:next_idx,cells_to_use(iC)));
    pks = findpeaks(zscore(ms.RawTraces(this_idx:next_idx,cells_to_use(iC))),33,'MinPeakDistance',2, 'MinPeakHeight', 5);
        d_max(iC) = max(this_data);
    d_pks(iC) = length(pks);
end
[best_pks, best_cells] = sort(d_pks, 'descend');
[best_most_pks, best_most_cells] = sort(d_pks(best_cells(1:40)), 'descend');
cells = best_cells(1:good_cells);
% cells = 1:good_cells;

sub_mat = reshape(1:((4+length(cells))*4),4,4+length(cells))'; % matrix to pull subplot values from
%%
if ishandle(108)
    close(108)
end
figure(108)
subplot(4+length(cells),4,1:16)
imagesc(b_f{1})
title(sprintf('Time = %.3fs', ms.time(this_idx)/1000));
set(gca, 'xtick', [], 'ytick', [])

set(gcf, 'color', 'k')
for iSub = 1:length(cells)
    this_sub = sub_mat(iSub+4,1:2);
    ax(iSub)= subplot(4+length(cells),4,sort(this_sub(:)));
    %
    this_data = zscore(ms.RawTraces(this_idx:next_idx,cells(iSub)));
    hold on
    plot(ms.time(this_idx:next_idx)/1000, this_data, 'color', c_ord(iSub,:), 'linewidth', 2);
    this_data(~ms.Binary(this_idx:next_idx,cells(iSub))) = NaN;
    h = area(ms.time(this_idx:next_idx)/1000,this_data,1);
    h.FaceColor = c_ord(iSub,:);
    h.EdgeColor = h.FaceColor;
    h.LineWidth = 0.01;
    h.BaseLine.Color =[0 0 0 0];
    % MS_Ca_Raster(ms.Binary(1:ms.timestamps(2),:)'); %, ms.time(1:ms.timestamps(2))'/1000
    ylim([-3 10])
    % ylabel(num2str(cells(iSub)), 'Rotation',0);
    if iSub == length(cells)
        set(gca, 'ytick', [], 'color', 'k');
        set(gca,'Color', 'none', 'YColor', 'none', 'XColor', 'w');
        ax1 = gca;
        ax1.XAxis.Label.Color=c_ord(1,:);
        ax1.XAxis.Label.FontWeight = 'Bold';
        ax1.XAxis.Label.Visible='on';
        
    else
        axis off
        set(gca, 'xtick', [], 'ytick', []);
    end
end
linkaxes(ax, 'x');

this_sub = sub_mat(5:end,3:4);
subplot(4+length(cells),4,sort(this_sub(:)))
imagesc(ca_f{1})
colormap('gray')
set(gca, 'xtick', [], 'ytick', [])
hold on
for iC = length(cells):-1:1
    [x(iC),y(iC)]=find(ms.SFPs(:,:,cells(iC)) == max(ms.SFPs(:,:,cells(iC)),[],'all'));
end

scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5)

%%
set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

clear F
for iF = Fs:4:size(b_t,2)-Fs
    
    subplot(4+length(cells),4,1:16)
    imagesc(b_f{iF})
    title(sprintf('Time = %.3fs', ms.time(this_idx+iF)/1000));
    set(gca, 'xtick', [], 'ytick', [])
    
    %     title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(iF), ms.frameNum(this_idx+iF)*Fs));
    this_sub = sub_mat(iSub+4,1:2);
    subplot(4+length(cells),4,sort(this_sub(:)))
    xlim([ms.time(this_idx+iF)/1000 - 5 ms.time(this_idx+iF)/1000])
    %     h = vline(median([ms.timestamps(1)-Fs + iF  ms.timestamps(1)+Fs + iF]));
    
    this_sub = sub_mat(5:end,3:4);
    subplot(4+length(cells),4,sort(this_sub(:)))
    imagesc(ca_f{iF})
    hold on
    scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5);
    set(gca, 'xtick', [], 'ytick', [])
    
    drawnow;
    F(iF) = getframe(gcf) ;
    delete(h);
    
end


% set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

%%
parts = strsplit(cd, filesep);
mkdir(['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\JC_inter\' parts{end-1}])
writerObj = VideoWriter(['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\JC_inter\' parts{end-1} filesep 'Example_Ca_behav_best.avi']);
writerObj.FrameRate = Fs/4;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = Fs:4:size(b_t,2)-Fs
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


