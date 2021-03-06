addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared')); 

inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2')); 
% addpath('C:\Users\ecarm\Documents\GitHub\OASIS_matlab')
% oasis_setup;
% disp('OASIS added for deconv');

% raw_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1069\10_18_2019_PV1069_HATD5\H10_M52_S0_REM86s';
raw_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LT&sleep\12_15_2019_540day5\H13_M51_S34rem'; 
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
ms_dir = raw_dir;
lfp_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\540\12_15_2019_540day5'; 
%%
clearvars -except raw_dir ms_dir lfp_dir
cd(raw_dir)

vid_num = 1;

% b_obj = VideoReader(['behavCam' num2str(vid_num) '.avi']);
% for iF = b_obj.NumFrames:-1:1
%     b_f{iF} = read(b_obj, iF);
%     b_t(iF) = b_obj.CurrentTime;
% end

ca_obj = VideoReader(['msCam' num2str(vid_num) '.avi']);
for iF = ca_obj.NumFrames:-1:1
    ca_f{iF} = read(ca_obj, iF);
    ca_t(iF) = ca_obj.CurrentTime;
end


%% get the ms struct
parts = strsplit(raw_dir, filesep); 

cd([lfp_dir filesep parts{end}])
fnames = dir('ms*.mat');
load(fnames.name)
ms = msExtractBinary_detrendTraces(ms_seg);
clear ms_seg
% frame_offset = ms.timestamps(1);
Fs = mode(diff(ms.time));

vid_idx = [1 find(diff(ms.frameNum) ~=1)+1]; % get the indicies where the video starts.
this_idx = vid_idx(vid_num);
if vid_num == length(vid_idx)
    next_idx = length(ms.RawTraces);
else
    next_idx = vid_idx(vid_num+1);
end

ms.time = ms.time - ms.time(1);
%% get the LFP data from the cut MS (if it is unchanged.  this is true for PV_1069_HATD5 H10_M52_S0_REM86s)

cd(lfp_dir);
load('ms_resize.mat');

this_seg = find(contains(ms_seg_resize.file_names, parts{end})); 
disp(ms_seg_resize.file_names{this_seg})

csc = ms_seg_resize.NLX_csc{this_seg};
if contains(ms_seg_resize.file_names{this_seg}, 'SW','IgnoreCase',true)
    csc_raw = csc;
    csc_raw.data = csc.data(2,:);
    csc_raw.label = csc.label{2};
    csc_raw.cfg.hdr = [];
    csc_raw.cfg.hdr{1} = csc.cfg.hdr{1};
    
    cfg_filt_d = [];
    cfg_filt_d.type = 'butter'; %Cheby1 is sharper than butter
    cfg_filt_d.f  = [140 250]; % broad, could use 150-200?
    cfg_filt_d.order = 4; %type filter order (fine for this f range)
    cfg_filt_d.display_filter = 0; % use this to see the fvtool
    ripple = FilterLFP(cfg_filt_d, csc_raw);
end

if contains(ms_seg_resize.file_names{this_seg}, 'REM','IgnoreCase',true)
    these_SWDs = ms_seg_resize.SWD_evts{this_seg};
    SWD_csc = csc.data(2,:);
    keep_idx = [];
    for iS = 1:length(these_SWDs.tstart)
        keep_idx = [keep_idx , (nearest_idx(these_SWDs.tstart(iS),csc.tvec) -1000)   : nearest_idx(these_SWDs.tend(iS),csc.tvec)+2000];
    end
    idx = zeros(1,length(SWD_csc));
    idx(unique(keep_idx)) = 1; 
    SWD_csc(~idx) = NaN;
end

csc.tvec = csc.tvec - csc.tvec(1); 
csc_idx = nearest_idx(ms.time(this_idx:next_idx)/1000,csc.tvec);

%% deconv debugger.
% % counter_init(size(ms.RawTraces,2),size(ms.RawTraces,2))
% for iChan = 25:-1:1
%     counter(iChan, size(ms.RawTraces,2))
%     [denoise,deconv] = deconvolveCa(ms.detrendRaw(:,iChan), 'foopsi', 'ar1', 'smin', -3, 'optimize_pars', true, 'optimize_b', true);
%     ms.denoise(:,iChan) = denoise;
%     ms.deconv(:,iChan) = deconv;
% end

% debugging
% if ishandle(111)
%     close(111)
% end
% figure(111)
% % iChan = 25;
% % hold on
% % plot(ms.detrendRaw(:,iChan), 'k');
% % plot(ms.denoise(:,iChan)-.2, 'r');
% % plot(((ms.deconv(:,iChan)./max(ms.deconv(:,iChan)))*.1) -.2, 'b');
% % plot((ms.Binary(:,iChan)*.1)+max(ms.RawTraces(:,iChan)), 'g');
% MS_plot_ca_trace(ms.FiltTraces(1:33*60,1:50)')

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
parts = strsplit(raw_dir, filesep); 
good_cells =20;
% find x nice cells.
if exist('SA.mat', 'file')
    load('SA.mat')
    cells_to_use = SA.WholePlaceCell; 
    c_ord = linspecer(length(cells));
% c_ord(((good_cells/2)-1):(good_cells/2)+4,:) = []; % remove yellows. They look terrible. 
else
    cells_to_use = 1:size(ms.RawTraces,2); 
    c_ord = linspecer(good_cells+6);
c_ord(((good_cells/2)-1):(good_cells/2)+4,:) = []; % remove yellows. They look terrible. 
end

d_max = []; 
for iC = length(cells_to_use):-1:1
       this_data = zscore(ms.RawTraces(this_idx:next_idx,cells_to_use(iC)));
       warning('OFF'); 
    pks = findpeaks(zscore(ms.RawTraces(this_idx:next_idx,cells_to_use(iC))),33,'MinPeakDistance',2, 'MinPeakHeight', 5);
        d_max(iC) = max(this_data);
    d_pks(iC) = length(pks);
end

[best_pks, best_cells] = sort(d_max, 'descend');
% [best_most_pks, best_most_cells] = sort(d_pks(best_cells(1:40)), 'descend');
% if exist(['C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\SeqNMF_EC\' parts{end-2} filesep parts{end-1} filesep 'all_REM_sweeps_parts.mat'], 'file')
%     load('all_REM_sweeps_parts.mat');
%     cells = all_sweeps.parts{1}.indSort(1:good_cells); 
% else
cells = 1:good_cells;
% end

sub_mat = reshape(1:((4+length(cells))*4),4,4+length(cells))'; % matrix to pull subplot values from
%%
if ishandle(108)
    close(108)
end
figure(108)
subplot(4+length(cells),4,1:16)
hold on
plot(csc.tvec(csc_idx(1):csc_idx(end)), csc.data(2,csc_idx(1):csc_idx(end)), 'color', c_ord(3,:), 'LineWidth', 1);
if contains(ms_seg_resize.file_names{this_seg}, 'REM','IgnoreCase',true)
    plot(csc.tvec(csc_idx(1):csc_idx(end)), SWD_csc(csc_idx(1):csc_idx(end)), 'color', c_ord(12,:), 'linewidth', 1);
%     plot(csc.tvec, csc.data(4,:), 'color', c_ord(7,:), 'LineWidth', 1); % plot theta filder
elseif contains(ms_seg_resize.file_names{this_seg}, 'SW','IgnoreCase',true)
        plot(ripple.tvec, ripple.data, 'color', c_ord(7,:), 'LineWidth', 1); % plot theta filder
end
y_val = [min(csc.data(2,:)) max(csc.data(2,:))];
xlim([csc.tvec(1) csc.tvec(csc_idx(end))])
% title(sprintf('Time = %.3fs', csc.tvec());
set(gca, 'xtick', [], 'ytick', [], 'color', 'k')

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
% xlim([0 700]); 
ylim([100 size(ca_f{1},1)])
set(gca, 'xtick', [], 'ytick', [])
hold on
for iC = length(cells):-1:1
    [x(iC),y(iC)]=find(ms.SFPs(:,:,cells(iC)) == max(ms.SFPs(:,:,cells(iC)),[],'all'));
end

scatter(y*ms.ds, x*ms.ds,125,c_ord,'LineWidth',1.5)

%%
set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

clear F
for iF = Fs:4:size(ca_t,2)-Fs
    
    subplot(4+length(cells),4,1:16)
    ylim(y_val); 
    xlim([csc.tvec(csc_idx(iF)) - 5 csc.tvec(csc_idx(iF))])
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
    ylim([100 size(ca_f{1},1)]);
    hold on
    scatter(y*ms.ds, x*ms.ds,125,c_ord,'LineWidth',1.5);
    set(gca, 'xtick', [], 'ytick', [])
    
    drawnow;
    F(iF) = getframe(gcf) ;
    delete(h);
    
end


% set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

%%
parts = strsplit(cd, filesep);
mkdir(['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\EV_inter\' parts{end-1}])
writerObj = VideoWriter(['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\EV_inter\' parts{end-1} filesep 'Example_Ca_' ms_seg_resize.file_names{this_seg} '.avi']);
writerObj.FrameRate = Fs/4;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = Fs:4:size(ca_t,2)-Fs
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


