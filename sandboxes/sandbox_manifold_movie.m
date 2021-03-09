%% manifold movie sandbox

addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2')); 

cd('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\HATD5')
inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold'; 
load('ms.mat');
load('v2_manifold.mat');
load('behav.mat');
clearvars -except ms behav v2_manifold

dir_parts = strsplit(cd,  filesep); 
session = dir_parts{end};
session = strrep(session, 'HATDS', 'HATS'); % in case of HATDSwitch vs HATSwitch.
subject = dir_parts{end-1};



cd(['C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\SeqNMF_EC\' subject])

this_dir = MS_list_dir_names(cd, session);
cd(this_dir{1})

Seq_REM = load('all_REM_sweeps_parts.mat', 'all_sweeps');
Seq_REM = Seq_REM.all_sweeps; 
Seq_wake = load('all_wake_seqs_parts.mat', 'all_seq');
Seq_wake = Seq_wake.all_seq; 
%% align behav is needed. 

[u_tvec, u_idx] = unique(behav.time);
pos(:,1) = interp1(u_tvec,behav.position(u_idx,1),ms.time);
pos(:,2) = interp1(u_tvec,behav.position(u_idx,2),ms.time);
velocity = interp1(u_tvec, behav.speed(u_idx), ms.time);
move_idx = velocity > 2.5; 
%% plot an example for track data.v

% actual data
this_type = v2_manifold.v2_track(move_idx,:); 
this_type = this_type(1:size(Seq_wake.parts.H,2),:);

% get seq indices. 
if sum(Seq_wake.parts.is_significant) ~=0
    fprintf('Sig seq detected in SeqNMF.  Using window = %is\n', Seq_wake.parts.L);
end
Seq_idx = Seq_wake.parts.H >= prctile(Seq_wake.parts.H, 99.5);

Seq_samples = Seq_wake.parts.L*mode(diff(ms.time)); % convert length of factor into samples. 

Seq_idx_out = zeros(size(Seq_idx)); 
for ii = length(Seq_idx):-1:1
    if Seq_idx(ii)
        Seq_idx_out(ii:ii+Seq_samples) = 1; 
    end
end

% get the SEq times
% Seq_per = prctile(Seq_wake.parts.H, 99); 

%or try only the top x H values
    [H_pks, H_pks_idx] = findpeaks(Seq_wake.parts.H, 'MinPeakDistance', Seq_samples);
    [H_top, top_idx] = sort(H_pks, 'descend');
    H_top_rel = H_top/max(H_top);
    H_pks_idx_sort = H_pks_idx(top_idx(1:5));

Seq_idx_out = zeros(size(Seq_wake.parts.H));
Color_idx = Seq_idx_out; 
for ii = 1:length(H_pks_idx_sort)
    Seq_idx_out(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = 1;
    Color_idx(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = ii;
end
Seq_idx_out = logical(Seq_idx_out); 
these_idx = find(Seq_idx_out == 1);

figure(1)
subplot(4,1,1)
% hold on
% plot(Seq_idx, 'r');
% plot(Seq_idx_out, 'b');
imagesc((0:length(Seq_wake.data)-1)/mode(diff(ms.time)),1:size(Seq_wake.data,1), Seq_wake.data(Seq_wake.parts.indSort,:).*([1:size(Seq_wake.data,1)]+(floor(size(Seq_wake.data,1)/2)))')
ylabel('cell ID (sorted)')

c_ord = parula(size(H_pks_idx_sort,2));
for ii = 1:length(H_pks_idx_sort)
    rectangle('position', [H_pks_idx_sort(ii)/mode(diff(ms.time)), -10,20,Seq_wake.parts.L],'facecolor', c_ord(ii,:), 'edgecolor', [1 1 1])
end
set(gca, 'color', c_ord(1,:))
ylim([-10 size(Seq_wake.data,1)])
% subplot(4,4,[5 9 13])
% title('Sig factor')
% ylabel('cell ID (sorted)')
% imagesc(squeeze(Seq_wake.parts.W(Seq_wake.parts.indSort,:)).*([1:size(Seq_wake.parts.W,1)]+(floor(size(Seq_wake.data,1)/2)))')

subplot(4,1,2:4)
h = scatter3(this_type(~Seq_idx_out,2), this_type(~Seq_idx_out,3), this_type(~Seq_idx_out,4), 'o'); 
h.MarkerEdgeColor = [.8 .8 .8];
h.MarkerEdgeAlpha = .2;
xlim([min(this_type(~Seq_idx_out,2)), max(this_type(~Seq_idx_out,2))]);
ylim([min(this_type(~Seq_idx_out,3)), max(this_type(~Seq_idx_out,3))]);
zlim([min(this_type(~Seq_idx_out,4)), max(this_type(~Seq_idx_out,4))]);

hold on
scatter3(this_type(Seq_idx_out,2), this_type(Seq_idx_out,3), this_type(Seq_idx_out,4), 14*ones(1,length(these_idx)),Color_idx(Seq_idx_out), 'filled')

xlim([min(this_type(~Seq_idx_out,2)), max(this_type(~Seq_idx_out,2))]);
ylim([min(this_type(~Seq_idx_out,3)), max(this_type(~Seq_idx_out,3))]);
zlim([min(this_type(~Seq_idx_out,4)), max(this_type(~Seq_idx_out,4))]);

    view(-39.394802263855850, 23.550427232944173);

    set(gcf, 'position', [600    50   896   890]);
    
mkdir([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session])
saveas(gcf, [inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'Seq_wake_manifold.fig'])
saveas(gcf,[inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'Seq_wake_manifold.png'])

% % H = get(gcf, 'children');
% H = findobj(gcf,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');
% 
% for iH = 1:length(H)
%     set(H(iH), 'fontsize', 18, 'fontname', 'Helvetica', 'TickDir', 'out', 'fontweight', 'normal')
% end

% this_dur = (33*60*2); % last value is number of minutes to run. 
% 
% F(this_dur) = getframe(gcf); % pre allocate; 

%% make a track data movie. 
pos_c_ord = parula(ceil(max(pos(:,1)))+1);
c_ord = parula(size(H_pks_idx_sort,2));
tvec = 1:length(this_type); 
tvec = tvec*(1/mode(diff(ms.time))); 
this_type_gpu = gpuArray(this_type); 

figure('Position',[680 485 531 493],'MenuBar','none','ToolBar','none','resize','off')
set(gcf, 'position', [440   185   560   613])
for ii = 4:(30*60*2)+4
    subplot(4,1,1:3)
        hold on
%     scatter3(v2(:,2), v2(:,3), v2(:,4), 6*ones(1,length(position_per_frame)),position_per_frame)
%     alpha(.2)
if Seq_idx_out(ii)
        c1 = plot3(this_type_gpu(ii-3:ii,2), this_type_gpu(ii-3:ii,3), this_type_gpu(ii-3:ii,4),'color', c_ord(Color_idx(ii),:), 'linewidth', 4);

    c3 = plot3(this_type_gpu(ii,2), this_type_gpu(ii,3), this_type_gpu(ii,4), '.', 'color',c_ord(Color_idx(ii),:), 'markersize', 10);
    
else
    c1 = plot3(this_type_gpu(ii-3:ii,2), this_type_gpu(ii-3:ii,3), this_type_gpu(ii-3:ii,4),'color', [0.8 .8 .8 .2], 'linewidth', 4);

    c3 = plot3(this_type_gpu(ii,2), this_type_gpu(ii,3), this_type_gpu(ii,4), 'o', 'color',[0.8 .8 .8 .2], 'markersize', 6);
    
end
    hold off
    xlim([min(this_type(:,2)) max(this_type(:,2))]);
    ylim([min(this_type(:,3)) max(this_type(:,3))]);
    zlim([min(this_type(:,4)) max(this_type(:,4))]);
    view(-39.394802263855850, 23.550427232944173);

    subplot(4,1,4)
    if Seq_idx_out(ii)
        plot(pos(ii,1),0,'.', 'color', pos_c_ord(round(pos(ii,1))+1,:), 'markersize', 30)
        hold on
        plot(pos(ii-2:ii,1),0, 'color', pos_c_ord(round(pos(ii,1))+1,:), 'markersize', 30, 'linewidth', 4)
        hold off
    else
        plot(pos(ii,1),0,'.', 'color', [.8 .8 .8 .2], 'markersize', 30)
        hold on
        plot(pos(ii-2:ii,1),0, 'color', [.8 .8 .8 .2], 'markersize', 30, 'linewidth', 4)
        hold off
    end
    % add time. 
    t1 = text(10, .5, [num2str(tvec(ii),3) 's']);
        
    xlim([min(pos(:,1)) max(pos(:,1))]);
    ylim([-1 1]); set(gca, 'ytick', []); 
    xlabel('position (cm)')
    set(gca, 'XDir','reverse')

     F(ii) = getframe(gcf) ;
      drawnow
      
      delete(c1);
      delete(t1); 
end

%% write to file. 
  writerObj = VideoWriter([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'wake_short.avi']);
  writerObj.FrameRate = round(33);
  writerObj.Quality = 80; 
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for ii=4:(30*60*2)+4
    % convert the image to a frame
    frame = F(ii) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% manifold movie with Seq using post_REM
% actual data
this_type = v2_manifold.v2_postREM(:); 

% get seq indices. 
for iSeq = 1:length(Seq_REM.parts)
    if sum(Seq_REM.parts{iSeq}.sig) ~=0
        fprintf('Sig seq detected in SeqNMF.  Using window = %is\n', Seq_REM.parts{iSeq}.L);
        good_L(iSeq) = iSeq;
    end
end

this_L = good_L(good_L ~=0);

Seq_idx = Seq_REM.parts{this_L}.H >= prctile(Seq_REM.parts{this_L}.H, 99.5);

Seq_samples = Seq_REM.parts{this_L}.L*mode(diff(ms.time)); % convert length of factor into samples. 

Seq_idx_out = zeros(size(Seq_idx)); 
for ii = length(Seq_idx):-1:1
    if Seq_idx(ii)
        Seq_idx_out(ii:ii+Seq_samples) = 1; 
    end
end

% get the SEq times
% Seq_per = prctile(Seq_wake.parts.H, 99); 

%or try only the top x H values
    [H_pks, H_pks_idx] = findpeaks(Seq_wake.parts.H, 'MinPeakDistance', Seq_samples);
    [H_top, top_idx] = sort(H_pks, 'descend');
    H_top_rel = H_top/max(H_top);
    H_pks_idx_sort = H_pks_idx(top_idx(1:5));

Seq_idx_out = zeros(size(Seq_wake.parts.H));
for ii = 1:length(H_pks_idx_sort)
    Seq_idx_out(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = 1;
end

these_idx = find(Seq_idx_out == 1);

figure(1)
subplot(4,1,1)
hold on
plot(Seq_idx, 'r');
plot(Seq_idx_out, 'b');


subplot(4,1,2:4)
h = scatter3(this_type(~Seq_idx_out,2), this_type(~Seq_idx_out,3), this_type(~Seq_idx_out,4), 'o'); 
h.MarkerEdgeColor = [.8 .8 .8];
h.MarkerEdgeAlpha = .2;

hold on
scatter3(this_type(Seq_idx_out,2), this_type(Seq_idx_out,3), this_type(Seq_idx_out,4), 6*ones(1,length(these_idx)),pos(Seq_idx_out,1), 'filled')


[az, el] = view;

this_dur = (33*60*2); % last value is number of minutes to run. 

F(this_dur) = getframe(gcf); % pre allocate; 





