function MS_JC_example_REM_video(raw_dir, ms_dir, decode_dir,output_dir, replay_idx)


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

raw_dir = 'K:\Jisoo_Project\RawData\pv1192\4_17_2021_PV1192_HATD1'; 
decode_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1192\HATD1'; 
ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1192\4_17_2021_PV1192_HATD1'; 
replay_idx = [2584 5030 5316 ]; 

output_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Decoding_data\Videos\EC_museum'

re_len = 16; 
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
    mouse_image = imread('C:\Users\ecarm\Pictures\mouse_sleep_fill.jpeg');
end
%% generate a time vector for the all_binary_post_REM

REM_idx = contains(ms_seg_resize.hypno_label, 'REM');
post_idx = contains(ms_seg_resize.pre_post, 'post');

REM_blocks = find(REM_idx & post_idx);

all_tvec = [];
all_ca_f = [];
ca_mat = [];
ca_t_all = [];
for iB = REM_blocks
    
    all_tvec = [all_tvec ; ms_seg_resize.time{iB}];
    
%     get the raw videos while you are here
    v_files = dir([raw_dir filesep ms_seg_resize.file_names{iB} filesep 'msCam*']);
    
    for iV = 1:length(v_files)
        ca_obj = [];
        ca_obj = VideoReader([v_files(iV).folder filesep v_files(iV).name]);
        for iF = ca_obj.NumFrames:-1:1
            ca_f{iF} = read(ca_obj, iF);
            ca_t(iF) = ca_obj.CurrentTime;
            ca_mat(:,:,iF) = im2double(ca_f{iF});
        end
        
        ca_t_all = [ca_t_all, ca_t];
        % try to substract mean per pixel
        %     ca_mat = cell2mat(ca_f);
        %     ca_mat = ca_mat - mean(ca_mat, 3);
        
        all_ca_f = cat(3,all_ca_f, ca_mat) ;
        clear ca_mat ca_f ca_t
    end
%         ca_mov{iB}.ca_f = ca_f;
%         ca_mov{iB}.ca_t = ca_t;
        ca_mov{iB}.ca_ms_t = ms_seg_resize.time{iB};
    
end

% all_ca_f(all_ca_f< 0) = 0; 

if size(all_tvec, 1) ~= size(all_binary_post_REM,1)
    error('ms_seg post_rem tvec is not the same size as the all_binary_post_REM')
end

% mean_f = mean(all_ca_f,3); 

% get the sampling freq
Fs = mode(diff(all_tvec));

%%  generate the plots
for iEvt = 1:length(replay_idx)
replay_evt = iEvt;
start_idx = replay_idx(replay_evt) - round(Fs*1);
end_idx = replay_idx(replay_evt) + round(Fs *1); 

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
% if contains(raw_dir, 'HAT')
%     rectangle('position', [0 0 track_size(2)/2 1], 'facecolor', [core_colors(2,:) .4], 'EdgeColor', [core_colors(3,:), 0]);
%     rectangle('position', [track_size(2)/2 0 track_size(2) 1], 'facecolor', [core_colors(3,:) .4], 'EdgeColor', [core_colors(3,:), 0])
% end
% xlim([0 100])
title(sprintf('Time = %.2fs', all_tvec(1)/1000));
set(gca, 'xtick', [], 'ytick', [])

% plot the decoder
% ax(100) = subplot(6,4,5:6);
% imagesc(ms.time/1000, 1:size(decoding.WAKE_decoded_probabilities,1), decoding.WAKE_decoded_probabilities)
% % caxis([0 max(decoding.REM_decoded_probabilities(3:end-3,:), [], 'all')])

set(gcf, 'color', 'k')

% ax(200)= subplot(6,4,[ 5 6 9 10 13 14 17 18 21 22]);
 ax(200)= subplot(7,4,[7 8 11 12 15 16 19 20 23 24 27 28]); 
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
% ax(300) = subplot(6,4,[7 8 11 12 ]);%15 16 19 20 23 24]); 
ax(300) = subplot(7,4,[ 5 6 9 10]);%15 16 19 20 23 24]); 
hold on
% yyaxis left
% set(gca, 'color', 'k', 'XColor',[0 0 0], 'YColor',[0 0 0] );
% yyaxis right
plot(all_tvec(start_idx-Fs*1.5:end_idx+Fs*1.5)/1000, decoding.REM_decoded_position(start_idx-Fs*1.5:end_idx+Fs*1.5), 'o', 'color', [.4 .4 .4 .3], 'MarkerSize', 2, 'MarkerFaceColor', [.4 .4 .4]);

set(gca, 'color', 'k', 'XColor',[.6 .6 .6], 'YColor',[.6 .6 .6] ); %set background color. 
ylabel('Decoded position (cm)', 'Color', [.6 .6 .6])
xlabel('time (s)', 'Color', [.6 .6 .6])
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


% plot bayes theorum? 
subplot(7,4,[13 17])
C = [0 .5] ;  % center 
a = .8;      % major axis 
e = .8 ;    % eccentricity 
b = a*sqrt(1-e^2) ; % minor axis 
th = linspace(0,2*pi) ; 
xe = C(1)+a*cos(th) ; 
ye = C(2)+b*sin(th) ; 
plot(xe,ye,'--w', 'LineWidth', 2);
text(0, .575, 'p (position | spikes)', 'color', 'w', 'fontsize', 10,'HorizontalAlignment', 'center')
hold on
% arrow
p1 = [0 1.2];                         % First Point
p2 = [1 2];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'color', 'w', 'LineWidth', 2, 'MaxHeadSize', 1)

% arrow ca 2 bay
p1 = [1.2 -1.2];                         % First Point
p2 = [0 -.2];                         % Second Point
dp = p2-p1;                         % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'color', 'w', 'LineWidth', 2, 'MaxHeadSize', 1)

% % flat
% p1 = [1.2 -1.2];                         % First Point
% p2 = [0 -1.2];                         % Second Point
% dp = p2-p1;                         % Difference
% quiver(p1(1),p1(2),dp(1),dp(2),1, 'color', 'w', 'LineWidth', 2, 'MaxHeadSize', 1, 'ShowArrowHead', 'off')

xlim([-.8 1.2])
ylim([-2 2.2])
    set(gca, 'color', 'k'); axis off; box on;



% add some calcium image
% subplot(6,4,[ 17 18 21 22]);
%     set(gca, 'color', 'k'); axis off; box on;
% hold on
% % text(
% xlim([0 size(track_image,2)])

% mouse sleeping

subplot(7,4,[ 21 22 25 26]);
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

for iF = start_idx:end_idx
    
%     htrack = subplot(6,4,1:4); 
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

           h.Marker = '<';
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
    subplot(7,4,[7 8 11 12 15 16 19 20 23 24 27 28]); 
    xlim([all_tvec(iF)/1000 - 1.5 all_tvec(iF)/1000 + 1.5])
        hl = xline(median([all_tvec(iF)/1000 - 1.5 all_tvec(iF)/1000 + 1.5]), 'linewidth', 6 );
    
%     subplot(6,4,[7 8 11 12 15 16 19 20 23 24])
%     imagesc((all_ca_f(:,:,iF) - mean_f)./mean_f)
% %     hold on
% %     scatter(y*ms.ds, x*ms.ds,100,c_ord,'LineWidth',1.5);
%     set(gca, 'xtick', [], 'ytick', [])

% move the decoding
%     subplot(6,4,[7 8 11 12 ]); %15 16 19 20 23 24])
    subplot(7,4,[ 5 6 9 10])
    xlim([all_tvec(iF)/1000 - 1.5 all_tvec(iF)/1000 + 1.5])
    h2 = xline(median([all_tvec(iF)/1000 - 1.5 all_tvec(iF)/1000 + 1.5]), 'linewidth', 6 );
    
    set(gca, 'xtick', [], 'ytick', [])
    
    drawnow;
    F(iF) = getframe(gcf) ;
%     delete(hh);
    if ~ismember(iF, replay_win(iEvt,:))%[replay_idx(1)+[0:5] replay_idx(2)+[0:5]]) || iF == replay_win(iEvt,1)
        set(findobj(gcf, 'type','Scatter'), 'Visible','off')
        delete(h);
    end
%     if iF == replay_win(iEvt,end) +1 || iF == replay_win(iEvt,1) -1
%         htrack = [];
%             htrack = subplot(6,4,1:4); 
%     imagesc(track_image)
%     end
    delete(hl)
    delete(h2)

    
end

%% write the video to disc. 
rate = 1; 

parts = strsplit([raw_dir filesep ms_seg_resize.file_names{iB}], filesep);
mkdir([output_dir filesep parts{end-1}])
writerObj = VideoWriter([output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_short_' num2str(rate) 'x_evt' num2str(replay_evt) '_fix.avi']);
fprintf('Saving video as: <strong>%s</strong>\n', [output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_short_' num2str(rate) 'x.avi'])
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
end



%% make a gif?

filename = [output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_short_' num2str(rate) 'x_evt' num2str(replay_evt) '_fix.gif'] 
mov1 = VideoReader([output_dir filesep parts{end-1} filesep 'Example_Ca_Post_REM_decode_raster_short_' num2str(rate) 'x_evt' num2str(replay_evt) '_fix.avi']);
vidFrames = read(mov1);
for n = 1:mov1.NumFrames
      [imind,cm] = rgb2ind(vidFrames(:,:,:,n),255);
      if n == 1;
          imwrite(imind,cm,filename,'gif','DelayTime',1/((Fs/rate))*2, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','DelayTime',1/((Fs/rate))*2,'WriteMode','append');
      end
end



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


