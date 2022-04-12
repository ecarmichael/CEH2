%% sandbox JC museum video. 

% addpath(genpath('/home/williamslab/Documents/Github/CEH2'));
% 
% 
% inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
% addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
% 
% pv1060 LTD5
% raw_dir = '/home/williamslab/Desktop/7_19_2019_PV1060_LTD5';
% decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/LTD5';
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
% ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/LTD5';
% replay_idx = [991 1001 1636 1646 2361]; 

% %pv1060 HATD5
% raw_dir = '/media/williamslab/Seagate Expansion Drive/Jisoo_Project/RawData/pv1060/11_26_2019_PV1060_HATSwitch';
% decode_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1060/HATDSwitch';
% % ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';
% ms_dir = '/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/HATDSwitch';
% replay_idx = [556 581 591 646 1206 1216 1486 1511 2291 2296 2306]; 

raw_dir = 'C:\Users\ecarm\Desktop'; 
decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1254\LTD1'; 
ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_13_2021_pv1254_LTD1'; 
replay_idx = [773 2516];

% raw_dir = 'C:\Users\ecarm\Desktop'; 
% decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\LTD1'; 
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1'; 
% replay_idx = [2766,3750,5847, 8171, 8336];

% raw_dir = 'C:\Users\ecarm\Desktop'; 
% decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1192\HATD1'; 
% ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1192\4_17_2021_PV1192_HATD1'; 
% replay_idx = [2584 5030 5316 ]; 

output_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Decoding_data\Videos\EC_museum'

re_len = 17; 
core_colors = linspecer(3); 


%% load everything

load([ms_dir filesep 'all_binary_post_REM'])
load([ms_dir filesep 'all_RawTraces_post_REM.mat'])
warning off
load([ms_dir filesep 'ms_resize.mat'])
warning on

load([decode_dir filesep 'decoding.mat'])

if isunix
track_image = imread(['/home/williamslab/Dropbox (Williams Lab)/Inter/pv1060/LTD5/Track_crop.jpg']);
else
    track_image = imread('C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\LTD5\Track_crop.jpg');
    mouse_image = imread('C:\Users\ecarm\Pictures\sleeping_mouse_back.jpeg');
%     mouse_image = imread('C:\Users\ecarm\Pictures\mouse_sleep_fill.jpeg');

end

GE_video = 'C:\Users\ecarm\Dropbox (Personal)\Jisoo and Eric\246_Barnes_Ca_only.avi'; 
GE_vid = VideoReader(GE_video);
GE_frames = flipud(read(GE_vid));
GE_frames = GE_frames(:,:,:,300:end); 


replay_evt =1;

Fs = mode(diff(ms_seg_resize.time{1}));

iEvt = replay_evt; 
start_idx = replay_idx(replay_evt) - round(Fs*1);
end_idx = replay_idx(replay_evt) + round(Fs *1); 
this_dir = 'L'; 
%%  get the tvec

REM_idx = contains(ms_seg_resize.hypno_label, 'REM');
post_idx = contains(ms_seg_resize.pre_post, 'post');

REM_blocks = find(REM_idx & post_idx);

all_tvec = [];

for iB = REM_blocks
    
    all_tvec = [all_tvec ; ms_seg_resize.time{iB}];
end
%%
if ishandle(108)
    close(108)
end
figure(108)
subplot(7,4,[1:4])
    set(gca, 'color', 'k'); axis off; box on;
hold on
imagesc(track_image)
xlim([0 size(track_image,2)])
track_size = get(gca, 'xlim'); 

title(sprintf('Time = %.2fs', all_tvec(1)/1000));
set(gca, 'xtick', [], 'ytick', [])


set(gcf, 'color', 'k')

% % ax(200)= subplot(6,4,[ 5 6 9 10 13 14 17 18 21 22]);
% %  ax(200)= subplot(7,4,[7 8 11 12 15 16  ]); 
ax(200) = subplot(7,4,[19 20 23 24 27 28]);
MS_Ca_Raster(all_binary_post_REM(start_idx-Fs*1.5:end_idx+Fs*1.5,:)', all_tvec(start_idx-Fs*1.5:end_idx+Fs*1.5)/1000, 5)
set(gca, 'color', 'k', 'XColor',[.6 .6 .6], 'YColor',[.6 .6 .6] ); %set background color. 
ylabel('Cell ID', 'Color', [.6 .6 .6])
xlabel('time (s)', 'Color', [.6 .6 .6])
if exist('replay_idx', 'var')
   hold on
   replay_win = [];
    for iR = replay_evt%1:length(replay_idx)
        replay_win(iR,:) = replay_idx(iR)+[0:re_len];
        if start_idx <  replay_idx(iR) < end_idx
        r = rectangle('position', [all_tvec(replay_idx(iR))/1000, size(all_binary_post_REM,2)+.5, re_len/Fs, 4], 'FaceColor', core_colors(1,:), 'EdgeColor', core_colors(1,:));
        ylim([0 size(all_binary_post_REM,2)+10.5])
        end        
    end
end

% 
 ax(300) =  subplot(7,4,[7 8 11 12 15 16]);
% ax(300) = subplot(4,6,[11 12 17 18 23 24]);%15 16 19 20 23 24]); 
% ax(300) = subplot(7,4,[ 5 6 9 10]);%15 16 19 20 23 24]); 
hold on
plot(all_tvec(start_idx-Fs*1.5:end_idx+Fs*1.5)/1000, decoding.REM_decoded_position(start_idx-Fs*1.5:end_idx+Fs*1.5), 'o', 'color', [.4 .4 .4 .3], 'MarkerSize', 2, 'MarkerFaceColor', [.4 .4 .4]);

set(gca, 'color', 'k', 'XColor',[.6 .6 .6], 'YColor',[.6 .6 .6], 'xtick', [] ); %set background color. 
ylabel('Decoded position (cm)', 'Color', [.6 .6 .6])
% xlabel('time (s)', 'Color', [.6 .6 .6])

if exist('replay_idx', 'var')
   hold on
   replay_win = [];
    for iR = replay_evt%1:length(replay_idx)
        replay_win(iR,:) = replay_idx(iR)+[0:re_len];
        if start_idx <  replay_idx(iR) < end_idx
            r_colors = winter(length(replay_win(iR,:))); 
            for ii = 1:length(replay_win(iR,:))
                plot(all_tvec(replay_win(iR,ii))/1000, decoding.REM_decoded_position(replay_win(iR,ii)), 'o', 'color', r_colors(ii,:), 'MarkerFaceColor', r_colors(ii,:));
            end

        r = rectangle('position', [all_tvec(replay_idx(iR))/1000, max(decoding.REM_decoded_position)+5, re_len/Fs, 1], 'FaceColor', core_colors(1,:), 'EdgeColor', core_colors(1,:));
        ylim([ min(decoding.REM_decoded_position) max(decoding.REM_decoded_position)+10])
        end        
    end
end
set(gca, 'YTick', [floor(min(decoding.REM_decoded_position)):10:ceil(max(decoding.REM_decoded_position))])
ylabel('Decoded position (cm)')
linkaxes(ax, 'x');

% 
% % plot bayes theorum? 
% subplot(7,4,[14, 18])
% C = [0 .5] ;  % center 
% a = .8;      % major axis 
% e = .8 ;    % eccentricity 
% b = a*sqrt(1-e^2) ; % minor axis 
% th = linspace(0,2*pi) ; 
% xe = C(1)+a*cos(th) ; 
% ye = C(2)+b*sin(th) ; 
% plot(xe,ye,'--w', 'LineWidth', 2);
% text(0, .575, 'p (position | spikes)', 'color', 'w', 'fontsize', 10,'HorizontalAlignment', 'center')
% hold on
% % arrow
% p1 = [0 1.2];                         % First Point
% p2 = [0 2];                         % Second Point
% dp = p2-p1;                         % Difference
% quiver(p1(1),p1(2),dp(1),dp(2),1,'color', 'w', 'LineWidth', 2, 'MaxHeadSize', 1)
% 
% % arrow ca 2 bay
% p1 = [1.2 -1.2];                         % First Point
% p2 = [0 -.2];                         % Second Point
% dp = p2-p1;                         % Difference
% quiver(p1(1),p1(2),dp(1),dp(2),1,'color', 'w', 'LineWidth', 2, 'MaxHeadSize', 1)
% 
% % % flat
% % p1 = [1.2 -1.2];                         % First Point
% % p2 = [0 -1.2];                         % Second Point
% % dp = p2-p1;                         % Difference
% % quiver(p1(1),p1(2),dp(1),dp(2),1, 'color', 'w', 'LineWidth', 2, 'MaxHeadSize', 1, 'ShowArrowHead', 'off')
% 
% xlim([-.8 1.2])
% ylim([-2 2.2])
%     set(gca, 'color', 'k'); axis off; box on;
% 
% 

% add some calcium image
% subplot(4,6,[ 8 9  14 15 ]);
subplot(7,4,[5 6 9 10 13 14]);
    set(gca, 'color', 'k'); axis off; box on;
    [imind,cm] = rgb2ind(GE_frames(:,:,:,1),255);
imshow(imind, cm); 

% 
% % add some calcium image
% subplot(7,4,[13 17]);
%     set(gca, 'color', 'k'); axis off; box on;
% imshow(G(1).cdata); 

% mouse sleeping
% subplot(4,6,[ 19 ]);
subplot(7,4,[17 18 21 22 25 26]); 
    set(gca, 'color', 'k'); axis off; box on;
hold on
imagesc(flipud(imcomplement(mouse_image)))
xlim([0 size(mouse_image,2)])


% imagesc(all_ca_f(:,:,start_idx))
% colormap('gray')
% set(gca, 'xtick', [], 'ytick', [])
% hold on
% for iC = length(cells):-1:1
%     [x(iC),y(iC)]=find(ms.SFPs(:,:,cells(iC)) == max(ms.SFPs(:,:,cells(iC)),[],'all'));
% end
% 
% scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5)

set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

%%
clear F
max_deco_prob = max(decoding.REM_decoded_probabilities,[], 'all');
max_pos = max(decoding.REM_decoded_position); 

indx = start_idx:0.5:end_idx; 

for iF = start_idx:end_idx
    
%     htrack = subplot(6,4,1:4); 
%     htrack = subplot(4,6,1:6); 
    htrack = subplot(7,4,1:4); 

%     imagesc(track_image)
    title(sprintf('Time = %.2fs', all_tvec(iF)/1000), 'Color', [.6 .6 .6]);
    set(gca, 'xtick', [], 'ytick', [])
    hold on
    if isnan(decoding.REM_decoded_position(iF))
         h = [];
%          h3 = [];
    else
        scale_pos = (decoding.REM_decoded_position(iF) - min(decoding.REM_decoded_position))/ (max(decoding.REM_decoded_position) - min(decoding.REM_decoded_position)); 
                h = scatter(scale_pos * track_size(2) , size(track_image,1)/2, 800, 'LineWidth', 2);

%         h = scatter((decoding.REM_decoded_position(iF)/max_pos)*(size(track_image,2)), size(track_image,1)/2, 800, 'LineWidth', 5);
%         h = scatter(decoding.REM_decoded_position(iF), 1, 800, 'LineWidth', 5);
        if ismember(iF, replay_win(iEvt,:)) %&& iF < 1001
%            h3 = plot([(scale_pos * track_size(2)) (scale_pos * track_size(2))+100 ], [size(track_image,1)/2 size(track_image,1)/2]);
           if strcmp(this_dir, 'L')
               h.Marker = '>';
           else
           h.Marker = '<';
           end
           h.MarkerEdgeColor = r_colors(find(replay_win(iEvt,:)== iF),:); 
           h.MarkerFaceColor = r_colors(find(replay_win(iEvt,:)== iF),:);% core_colors(3,:);
%         elseif ismember(iF, replay_win(iEvt,:)) %&& iF >1000
%             h.MarkerEdgeColor = core_colors(3,:);
        else
           h.MarkerEdgeColor = [.4 .4 .4];
           h.Marker = '.';
        end
%         h.MarkerEdgeAlpha  =max(decoding.WAKE_decoded_probabilities(:,iF))/max_deco_prob;
    end
    
    % move the Calcium
    %     title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(iF), ms.frameNum(this_idx+iF)*Fs));
%     subplot(6,4,[5 6 9 10 13 14 17 18 21 22])
%     subplot(7,4,[7 8 11 12 15 16 19 20 23 24 27 28]);
% subplot(4,6,[11 12 17 18 23 24]);%15 16 19 20 23 24]); 
% subplot(7,4,[7 8 11 12 15 16])
subplot(7,4,[7 8 11 12 15 16])
    xlim([all_tvec(iF)/1000 - 1.5 all_tvec(iF)/1000 + 1.5])
        hl = xline(median([all_tvec(iF)/1000 - 1.5 all_tvec(iF)/1000 + 1.5]), 'linewidth', 6 );
    
%     subplot(6,4,[7 8 11 12 15 16 19 20 23 24])
%     imagesc((all_ca_f(:,:,iF) - mean_f)./mean_f)
% %     hold on
% %     scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5);
%     set(gca, 'xtick', [], 'ytick', [])

% move the decoding
%     subplot(6,4,[7 8 11 12 ]); %15 16 19 20 23 24])
%     subplot(7,4,[ 5 6 9 10])

%     if iF == replay_win(iEvt,end) +1 || iF == replay_win(iEvt,1) -1
%         htrack = [];
%             htrack = subplot(6,4,1:4); 
%     imagesc(track_image)
%     end
 

    % add some calcium image
% subplot(4,6,[8 9 10 13 14 15 16]);
subplot(7,4,[5 6 9 10 13 14]);
    set(gca, 'color', 'k'); axis off; box on;
    [imind,cm] = rgb2ind(GE_frames(:,:,:,find(indx == iF)+920),255);
imshow(imind, cm); 
    
 drawnow;
    F(iF) = getframe(gcf) ;
    delete(hl)
       if ~ismember(iF, replay_win(iEvt,:))%[replay_idx(1)+[0:5] replay_idx(2)+[0:5]]) || iF == replay_win(iEvt,1)
        set(findobj(gcf, 'type','Scatter'), 'Visible','off')
        delete(h);
    end
end

%% write the video to disc. 
rate = 1; 

parts = strsplit([ms_dir filesep ms_seg_resize.file_names{iB}], filesep);
mkdir([output_dir filesep parts{end-1}])
writerObj = VideoWriter([output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_museum_' parts{end-1} '_'  num2str(replay_evt) '.avi']);
fprintf('Saving video as: <strong>%s</strong>\n', [output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_museum_' parts{end-1} '_'  num2str(replay_evt) '.avi'])
writerObj.FrameRate = ((Fs)/rate)/4;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = start_idx:rate:end_idx
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);



%% make a gif?

filename = [output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_museum_' parts{end-1} '_'  num2str(replay_evt) '.avi'];
mov1 = VideoReader(filename);
vidFrames = read(mov1);
for n = 1:mov1.NumFrames
      [imind,cm] = rgb2ind(vidFrames(:,:,:,n),255);
      if n == 1
          imwrite(imind,cm,[filename(1:end-4) '.gif'],'gif','DelayTime',1/((Fs/rate))*5, 'Loopcount',inf);
      else
          imwrite(imind,cm,[filename(1:end-4) '.gif'],'gif','DelayTime',1/((Fs/rate))*5,'WriteMode','append');
      end
end



%% gif of the color cell activity

for n = 1:length(E)%ms.vidObj{1}.NumFrames-50
          [imind,cm] = rgb2ind(E(n).cdata,255);

      if n == 1
          imwrite(imind,cm,'temp_col.gif','gif','DelayTime',1/Fs, 'Loopcount',inf);
      else
          imwrite(imind,cm,'temp_col.gif','gif','DelayTime',1/Fs,'WriteMode','append');
      end
end


%% make a visualization of decoding
c_ord = winter(10); 
pos = 0:1:100; 
field = zeros(size(pos)); 

field([1 2 6 7]) = 1; 
field([3 5]) = 2; 
field([4]) = 3; 

figure(101)
clf
hold on
for ii = 1:10
sm_field(ii,:) = imgaussfilt(circshift(field, ii*9), round(MS_randn_range(1, 1, 2, 10)));
end

for ii = 1:10
h(ii) = area(pos, sm_field(ii,:)./max(sm_field,[], 'all'), 'FaceColor', c_ord(ii,:), 'FaceAlpha', .3)
end
legend([h(1), h(4) h(8)], {'cell 1', 'cell 2', 'cell 3'}, 'Box', 'off')
xlabel('position on track')
set(gca, 'ytick', [])

SetFigure([], gcf)

%%
% raster
figure(102)
clf
hold on
for ii = 1:10
%     this_data = smooth(all_RawTraces_post_REM(:,ii), 'sgolay'); 
   
% plot(all_tvec, (this_data./max(this_data))+ii)
this_s = find(sm_field(ii,:) >.25);
this_s = datasample(this_s, floor(length(this_s)/4)); 
line([this_s; this_s], [zeros(size(this_s))+ii-.5 ; zeros(size(this_s))+ii+.5], 'color', c_ord(ii,:), 'linewidth', 2)
    
end
% ylim([.5 10.5]*10)
set(gca, 'xticklabels', get(gca, 'xtick')/10);
set(gca, 'ytick', 1:10);
ylim([0.5 10.5])
xlabel('time (s)');
ylabel('neuron');
SetFigure([], gcf)
%%
% cd(raw_dir)
% 
% vid_num = 3;
% 
% b_obj = VideoReader(['behavCam' num2str(vid_num) '.avi']);
% for iF = b_obj.NumFrames:-1:1
%     b_f{iF} = read(b_obj, iF);
%     b_t(iF) = b_obj.CurrentTime;
% end
% 
% ca_obj = VideoReader(['msCam' num2str(vid_num) '.avi']);
% for iF = ca_obj.NumFrames:-1:1
%     ca_f{iF} = read(ca_obj, iF);
%     ca_t(iF) = ca_obj.CurrentTime;
% end
% 
% %%
% cd(ms_dir)
% load('ms.mat')
% ms = msExtractBinary_detrendTraces(ms);
% % frame_offset = ms.timestamps(1);
% Fs = mode(diff(ms.time));
% 
% vid_idx = [1 find(diff(ms.frameNum) ~=1)+1]; % get the indicies where the video starts.
% this_idx = vid_idx(vid_num);
% if vid_num == length(vid_idx)
%     next_idx = length(ms.RawTraces);
% else
%     next_idx = vid_idx(vid_num+1);
% end
% %% make the plot
% if ishandle(108)
%     close(108)
% end
% figure(108)
% subplot(6,4,[1:4])
%     set(gca, 'color', 'k'); axis off; box on;
% hold on
% imagesc(b_f{1})
% if contains(raw_dir, 'HAT')
%     rectangle('position', [0 0 50 1], 'facecolor', [core_colors(2,:) .4], 'EdgeColor', [core_colors(3,:), 0]);
%     rectangle('position', [50 0 50 1], 'facecolor', [core_colors(3,:) .4], 'EdgeColor', [core_colors(3,:), 0])
% end
% title(sprintf('Time = %.2fs', ms.time(1)/1000));
% set(gca, 'xtick', [], 'ytick', [])
% 
% % plot the decoder
% % ax(100) = subplot(6,4,5:6);
% % imagesc(ms.time/1000, 1:size(decoding.WAKE_decoded_probabilities,1), decoding.WAKE_decoded_probabilities)
% % % caxis([0 max(decoding.REM_decoded_probabilities(3:end-3,:), [], 'all')])
% 
% set(gcf, 'color', 'k')
% 
% ax(200)= subplot(6,4,[ 5 6 9 10 13 14 17 18 21 22]);
% MS_Ca_Raster(ms.Binary(this_idx:next_idx,:)', ms.time(this_idx:next_idx)/1000)
% linkaxes(ax, 'x');
% set(gca, 'color', 'k', 'XColor',[.6 .6 .6], 'YColor',[.6 .6 .6] ); %set background color. 
% ylabel('Cell ID', 'Color', [.6 .6 .6])
% xlabel('time (ms)', 'Color', [.6 .6 .6])
% 
% subplot(6,4,[7 8 11 12 15 16 19 20 23 24])
% imagesc(ca_f{1})
% colormap('gray')
% set(gca, 'xtick', [], 'ytick', [])
% hold on
% % for iC = length(cells):-1:1
% %     [x(iC),y(iC)]=find(ms.SFPs(:,:,cells(iC)) == max(ms.SFPs(:,:,cells(iC)),[],'all'));
% % end
% % 
% % scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5)
% 
% %%
% set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')
% max_pos = max(decoding.WAKE_decoded_position); 
% max_deco_prob = max(decoding.WAKE_decoded_probabilities,[], 'all');
% 
% clear F
% for iF = Fs:4:size(b_t,2)-Fs
%     
%     subplot(6,4,1:4)
%     hh = imagesc(b_f{iF});
%     title(sprintf('Time = %.3fs', ms.time(this_idx+iF)/1000));
%     set(gca, 'xtick', [], 'ytick', [])
%     hold on
%     if isnan(decoding.WAKE_decoded_position(iF))
%          h = [];
%     else
%         h = scatter((decoding.WAKE_decoded_position(this_idx+iF)/max_pos)*(size(b_f{iF},2)), size(b_f{iF},1)/2, 800, 'LineWidth', 5);
%         h.MarkerEdgeColor = core_colors(2,:);
% %         h.MarkerEdgeAlpha  =max(decoding.WAKE_decoded_probabilities(:,iF))/max_deco_prob;
%     end
%     %     title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(iF), ms.frameNum(this_idx+iF)*Fs));
%     subplot(6,4,[5 6 9 10 13 14 17 18 21 22])
%     xlim([ms.time(this_idx+iF)/1000 - 2.5 ms.time(this_idx+iF)/1000 + 2.5])
%         hl = vline(median([ms.time(this_idx+iF)/1000 - 2.5 ms.time(this_idx+iF)/1000 + 2.5]));
%     
%     subplot(6,4,[7 8 11 12 15 16 19 20 23 24])
%     imagesc(ca_f{iF})
% %     hold on
% %     scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5);
%     set(gca, 'xtick', [], 'ytick', [])
%     
%     drawnow;
%     F(iF) = getframe(gcf) ;
%     delete(hh);
%     delete(h); 
%     delete(hl)
%     
% end
% 

%% softmax the decoded probs
% rng(11)
% for ii = size(decoding.WAKE_decoded_probabilities,2):-1:1
%     if isnan(decoding.WAKE_decoded_probabilities(:,ii))
%         soft_decode_pos(ii) = NaN;
%     else
%         
%     soft_decode_pos(ii) = max(softmax(decoding.WAKE_decoded_probabilities(1:end-1,ii))); 
%     end
% end


