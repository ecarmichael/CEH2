%% manifold movie sandbox

inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

cd('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\HATD5')
% cd('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\LTD5')

inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
load('ms.mat');
load('v2_manifold.mat');
load('behav.mat');
clearvars -except ms behav v2_manifold inter_dir

dir_parts = strsplit(cd,  filesep);
session = dir_parts{end};
session = strrep(session, 'HATDS', 'HATS'); % in case of HATDSwitch vs HATSwitch.
subject = dir_parts{end-1};

% get the LFP
cd(['C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\' subject])
this_dir = MS_list_dir_names(cd, session);
cd([this_dir{1} filesep 'LFP_mats'])
REM_LFP = load('all_t_post_REM.mat', 'all_t_post_REM');
REM_LFP = REM_LFP.all_t_post_REM;


% get the seq data
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
data_to_plot = [nan(20,size(Seq_wake.data,2)); Seq_wake.data]; % just for plotting.  fills in with some nans
tvec = ms.time(move_idx)/1000; 
tvec = tvec(1:size(Seq_wake.parts.H,2)); % corrected tvec, which is odd because we have excluded imobility here.


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
subplot(4,2,1:2)
% hold on
% plot(Seq_idx, 'r');
% plot(Seq_idx_out, 'b');
hold on
% imagesc(data_to_plot)
imagesc((0:length(Seq_wake.data)-1)/mode(diff(ms.time)),1:size(Seq_wake.data,1), Seq_wake.data(Seq_wake.parts.indSort,:).*([1:size(Seq_wake.data,1)]+(floor(size(Seq_wake.data,1)/2)))')
ylabel('cell ID (sorted)')
set(gca,'YDir','reverse')
% axis xy
xlim([0 (length(Seq_wake.data)-1)/mode(diff(ms.time))]);
% box off
c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
for ii = 1:length(H_pks_idx_sort)
    rectangle('position', [H_pks_idx_sort(ii)/mode(diff(ms.time)), -20,Seq_wake.parts.L,20],'facecolor', c_ord_l(ii,:), 'edgecolor', [1 1 1])
end
% set(gca, 'color', c_ord(1,:))
ylim([-20 size(Seq_wake.data,1)])
% subplot(4,4,[5 9 13])
% title('Sig factor')
% ylabel('cell ID (sorted)')
% imagesc(squeeze(Seq_wake.parts.W(Seq_wake.parts.indSort,:)).*([1:size(Seq_wake.parts.W,1)]+(floor(size(Seq_wake.data,1)/2)))')

subplot(4,2,[3 5 7])
h = scatter3(this_type(~Seq_idx_out,2), this_type(~Seq_idx_out,3), this_type(~Seq_idx_out,4), 'o');
h.MarkerEdgeColor = [.8 .8 .8];
h.MarkerEdgeAlpha = .2;

hold on
% scatter3(this_type(Seq_idx_out,2), this_type(Seq_idx_out,3), this_type(Seq_idx_out,4), 14*ones(1,length(these_idx)),Color_idx(Seq_idx_out), 'filled')
for ii = 1:max(Color_idx)
    these_idx = Seq_idx_out & Color_idx == ii;
    h = scatter3(this_type(these_idx,2), this_type(these_idx,3), this_type(these_idx,4), 14*ones(1,length(this_type(these_idx,4))),repmat(c_ord_l(ii,:),size(this_type(these_idx,4))), 'filled');
end
legend({'non-seq', 'Seq 1', 'Seq 2', 'Seq 3', 'Seq 4', 'Seq 5'}, 'location', 'north', 'orientation', 'horizontal');
xlim([min(this_type(:,2)), max(this_type(:,2))]);
ylim([min(this_type(:,3)), max(this_type(:,3))]);
zlim([min(this_type(:,4)), max(this_type(:,4))]);
grid off

subplot(4,2,[4 6 8 ])
hold on
% scatter3(this_type(Seq_idx_out,2), this_type(Seq_idx_out,3), this_type(Seq_idx_out,4), 14*ones(1,length(these_idx)),Color_idx(Seq_idx_out), 'filled')
for ii = 1:max(Color_idx)
    these_idx = Seq_idx_out & Color_idx == ii;
    h = scatter3(this_type(these_idx,2), this_type(these_idx,3), this_type(these_idx,4), 12*ones(1,length(this_type(these_idx,4))),repmat(c_ord_l(ii,:),size(this_type(these_idx,4))), 'filled');
end
xlim([min(this_type(:,2)), max(this_type(:,2))]);
ylim([min(this_type(:,3)), max(this_type(:,3))]);
zlim([min(this_type(:,4)), max(this_type(:,4))]);


view(-39.394802263855850, 23.550427232944173);

set(gcf, 'position', [128 56 1446 890]);

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
% tvec = 1:length(this_type);
% tvec = tvec*(1/mode(diff(ms.time)));

tic
if ishandle(101)
    close(101)
end
figure(101)

subplot(7,1,1:4)
    xlim([min(this_type(:,2)) max(this_type(:,2))]);
    ylim([min(this_type(:,3)) max(this_type(:,3))]);
    zlim([min(this_type(:,4)) max(this_type(:,4))]);
    view(-39.394802263855850, 23.550427232944173);

subplot(7,1,6:7)
imagesc(1:size(Seq_wake.data,2),1:size(Seq_wake.data,1), Seq_wake.data(Seq_wake.parts.indSort,:).*([1:size(Seq_wake.data,1)]+(floor(size(Seq_wake.data,1)/2)))')
ylabel('cell ID (sorted)')
set(gca,'YDir','reverse')
% axis xy
xlim([min(tvec) max(tvec)]);
axis off
% xlabel('time (s)')

set(gcf, 'position', [440   185   560   660],'MenuBar','none','ToolBar','none','resize','off')
clear F
% for ii = 4:2:(33*60*.5)+4
for ii =3600:2:5282

    subplot(7,1,1:4)
    hold on
    %     scatter3(v2(:,2), v2(:,3), v2(:,4), 6*ones(1,length(position_per_frame)),position_per_frame)
    %     alpha(.2)
    if Seq_idx_out(ii)
        c1 = plot3(this_type(ii-3:ii,2), this_type(ii-3:ii,3), this_type(ii-3:ii,4),'color', c_ord_l(Color_idx(ii),:), 'linewidth', 4);
        
        c3 = plot3(this_type(ii,2), this_type(ii,3), this_type(ii,4), '.', 'color',c_ord_l(Color_idx(ii),:), 'markersize', 10);
        
    else
        c1 = plot3(this_type(ii-3:ii,2), this_type(ii-3:ii,3), this_type(ii-3:ii,4),'color', [0.8 .8 .8 .4], 'linewidth', 4);
        
        c3 = plot3(this_type(ii,2), this_type(ii,3), this_type(ii,4), 'o', 'color',[0.8 .8 .8 .2], 'markersize', 4);
        
    end
    hold off

    
    subplot(7,1,5)
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
    
    
    subplot(7,1,6:7)
%     if ii< 34
%         xlim([0 66])
%     else
        xlim([ii ii+66])
%     end
    
    F(ii) = getframe(gcf) ;
    drawnow
    
    delete(c1);
    delete(t1);
end
toc

%% write to file.
writerObj = VideoWriter([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'wake_short.avi']);
writerObj.FrameRate = 16.5;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
% for ii=4:2:(30*60*.2)+4
for ii =3962:2:5282

    % convert the image to a frame
    frame = F(ii) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% manifold movie with Seq using post_REM
% actual data

% get seq indices.
for iSeq = 1:length(Seq_REM.parts)
    if sum(Seq_REM.parts{iSeq}.sig) ~=0
        fprintf('Sig seq detected in SeqNMF.  Using window = %is\n', Seq_REM.parts{iSeq}.L);
        good_L(iSeq) = iSeq;
    end
end

this_L = good_L(good_L ~=0);


% use the sig L value
this_type = v2_manifold.v2_postREM;
this_type = this_type(1:length(Seq_REM.parts{this_L}.Train),:);

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
[H_pks, H_pks_idx] = findpeaks(Seq_REM.parts{this_L}.H, 'MinPeakDistance', Seq_samples);
[H_top, top_idx] = sort(H_pks, 'descend');
H_top_rel = H_top/max(H_top);
H_pks_idx_sort = H_pks_idx(top_idx(1:5));

Seq_idx_out = zeros(size(Seq_REM.parts{this_L}.H));
Color_idx = Seq_idx_out;
for ii = 1:length(H_pks_idx_sort)
    Seq_idx_out(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = 1;
    Color_idx(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = ii;
end
Seq_idx_out = logical(Seq_idx_out);
these_idx = find(Seq_idx_out == 1);

figure(2)
subplot(4,2,1:2)
% hold on
% plot(Seq_idx, 'r');
% plot(Seq_idx_out, 'b');
imagesc((0:length(Seq_REM.parts{this_L}.Train)-1)/mode(diff(ms.time)),1:size(Seq_REM.parts{this_L}.Train,1), Seq_REM.parts{this_L}.Train(Seq_REM.parts{this_L}.indSort,:).*([1:size(Seq_REM.parts{this_L}.Train,1)]+(floor(size(Seq_REM.parts{this_L}.Train,1)/2)))')
ylabel('cell ID (sorted)')
hv = vline(Seq_REM.parts{this_L}.seg_idx/mode(diff(ms.time)));
for ii = 1:length(hv)
    hv(ii).LineWidth = 2;
    hv(ii).Color = [.9 .9 .9];
end

c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
for ii = 1:length(H_pks_idx_sort)
    rectangle('position', [H_pks_idx_sort(ii)/mode(diff(ms.time)), -20,Seq_REM.parts{this_L}.L,20],'facecolor', c_ord_l(ii,:), 'edgecolor', [1 1 1])
end
set(gca, 'color', c_ord(1,:))
ylim([-20 size(Seq_REM.parts{this_L}.Train,1)])
box off
% subplot(4,4,[5 9 13])
% title('Sig factor')
% ylabel('cell ID (sorted)')
% imagesc(squeeze(Seq_wake.parts.W(Seq_wake.parts.indSort,:)).*([1:size(Seq_wake.parts.W,1)]+(floor(size(Seq_wake.data,1)/2)))')

subplot(4,2,[3 5 7])
h = scatter3(this_type(~Seq_idx_out,2), this_type(~Seq_idx_out,3), this_type(~Seq_idx_out,4), 'o');
h.MarkerEdgeColor = [.8 .8 .8];
h.MarkerEdgeAlpha = .2;

hold on
for ii = 1:max(Color_idx)
    these_idx = Seq_idx_out & Color_idx == ii;
    h = scatter3(this_type(these_idx,2), this_type(these_idx,3), this_type(these_idx,4), 14*ones(1,length(this_type(these_idx,4))),repmat(c_ord_l(ii,:),size(this_type(these_idx,4))), 'filled');
end
legend({'non-seq', 'Seq 1', 'Seq 2', 'Seq 3', 'Seq 4', 'Seq 5'}, 'location', 'north', 'orientation', 'horizontal');
xlim([min(this_type(:,2)), max(this_type(:,2))]);
ylim([min(this_type(:,3)), max(this_type(:,3))]);
zlim([min(this_type(:,4)), max(this_type(:,4))]);
view(-39.394802263855850, 23.550427232944173);
grid off


subplot(4,2,[4 6 8])
hold on
for ii = 1:max(Color_idx)
    these_idx = Seq_idx_out & Color_idx == ii;
    h = scatter3(this_type(these_idx,2), this_type(these_idx,3), this_type(these_idx,4), 14*ones(1,length(this_type(these_idx,4))),repmat(c_ord_l(ii,:),size(this_type(these_idx,4))), 'filled');
end
xlim([min(this_type(:,2)), max(this_type(:,2))]);
ylim([min(this_type(:,3)), max(this_type(:,3))]);
zlim([min(this_type(:,4)), max(this_type(:,4))]);
view(-39.394802263855850, 23.550427232944173);

set(gcf, 'position', [128 56 1446 890]);

% mkdir([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session])
saveas(gcf, [inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'Seq_postREM_manifold.fig'])
saveas(gcf,[inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'Seq_postREM_manifold.png'])


%% REM movie
this_lfp = REM_LFP(1:length(this_type)); % limit to training set from Seq data.
this_lfp = this_lfp./max(this_lfp);

c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
tvec = 1:length(this_type);
tvec = tvec*(1/mode(diff(ms.time)));

this_data = Seq_REM.parts{this_L}.Train; 

if ishandle(202)
    close(202)
end
figure(202)
set(gcf, 'position', [440   185   560   613],'MenuBar','none','ToolBar','none','resize','off')
subplot(8,1,1:5)
    xlim([min(this_type(:,2)) max(this_type(:,2))]);
    ylim([min(this_type(:,3)) max(this_type(:,3))]);
    zlim([min(this_type(:,4)) max(this_type(:,4))]);
    view(-39.394802263855850, 23.550427232944173);


subplot(8,1, 7:8)
imagesc((0:length(this_data)-1)/mode(diff(ms.time)),1:size(this_data,1), this_data(Seq_REM.parts{this_L}.indSort,:).*([1:size(this_data,1)]+(floor(size(this_data,1)/2)))')
ylabel('cell ID (sorted)')
set(gca,'YDir','reverse')
% axis xy
xlim([0 (length(this_data)-1)/mode(diff(ms.time))]);
% xlabel('time (s)')

c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
for ii = 1:length(H_pks_idx_sort)
    rectangle('position', [H_pks_idx_sort(ii)/mode(diff(ms.time)), -20,Seq_REM.parts{this_L}.L,20],'facecolor', c_ord_l(ii,:), 'edgecolor', [1 1 1])
end
set(gca, 'color', c_ord(1,:))
ylim([-20 size(Seq_REM.parts{this_L}.Train,1)])
box off



subplot(8,1,6)
imagesc(tvec, 1, this_lfp)
axis off
y_lim = ylim;
xlabel('time (s)')


clear F
for ii = 4:2:length(Seq_idx_out)
% for ii = 4:2:(30*60*2)+4
    subplot(8,1,1:5)
    hold on
    %     scatter3(v2(:,2), v2(:,3), v2(:,4), 6*ones(1,length(position_per_frame)),position_per_frame)
    %     alpha(.2)
    if Seq_idx_out(ii)
        c1 = plot3(this_type(ii-3:ii,2)', this_type(ii-3:ii,3)', this_type(ii-3:ii,4)','color', c_ord_l(Color_idx(ii),:)', 'linewidth', 4);
        
        c3 = plot3(this_type(ii,2)', this_type(ii,3)', this_type(ii,4)', '.', 'color',c_ord_l(Color_idx(ii),:)', 'markersize', 10);
        
    else
        c1 = plot3(this_type(ii-3:ii,2)', this_type(ii-3:ii,3)', this_type(ii-3:ii,4)','color', [0.8 .8 .8 .4], 'linewidth', 4);
        
        c3 = plot3(this_type(ii,2)', this_type(ii,3)', this_type(ii,4)', '.', 'color',[0.8 .8 .8 .2], 'markersize', 4);
        
    end
    
    hold off

    
    subplot(8,1,6)
    if ii< 34
        xlim([tvec(1) tvec(66)])
    else
        xlim([tvec(ii)-1 tvec(ii)+1])
    end
    ylim(y_lim);
    % add time.
    %     t1 = text(10, .5, [num2str(tvec(ii),3) 's']);
    %
    %     xlim([min(pos(:,1)) max(pos(:,1))]);
    %     ylim([-1 1]); set(gca, 'ytick', []);
    %     xlabel('position (cm)')
    %     set(gca, 'XDir','reverse')
    subplot(8,1,7:8)
%     if ii< 34
%         xlim([tvec(1) tvec(66)])
%     else
        xlim([tvec(ii) tvec(ii)+2])
%     end
    
    F(ii) = getframe(gcf) ;
    drawnow
    
    delete(c1);
    delete(t1);
end


%% write to file.
writerObj = VideoWriter([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'REM_full.avi']);
writerObj.FrameRate = round(16);
writerObj.Quality = 90;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
% for ii=4:2:(30*60*2)+4
for ii = 4:2:length(F)
    % convert the image to a frame
    frame = F(ii) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


