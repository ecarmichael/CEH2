%% manifold movie sandbox
addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2')); 


inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

% cd('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\HATD5')
% cd('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1069\LTD5')

% eva data
cd('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\Eva''s data\540\d1')

inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
load('ms_trk.mat'); % for eva data
ms_trk.time = ms_trk.time - ms_trk.time(1);
ms = ms_trk;
load('v2_manifold.mat');
load('behav.mat');
clearvars -except ms behav v2_manifold inter_dir

dir_parts = strsplit(cd,  filesep);
session = dir_parts{end};
session = strrep(session, 'd', 'day');
subject = dir_parts{end-1};

% get the SWD and SWR data
pro_data_dir = ['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\final_analysis' filesep subject];
this_dir = MS_list_dir_names(pro_data_dir, session);
cd([pro_data_dir filesep this_dir{1}])

SWD = load('swd_2.mat'); 
ripple = load('ripple_2.mat');


% get the LFP
% cd(['C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\' subject])
% this_dir = MS_list_dir_names(cd, session);
% cd([this_dir{1} filesep 'LFP_mats'])
% REM_LFP = load('all_t_post_REM.mat', 'all_t_post_REM');
% REM_LFP = REM_LFP.all_t_post_REM;


% get the seq data
cd(['C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\SeqNMF_EC\' subject])

this_dir = MS_list_dir_names(cd, session);
cd(this_dir{1})

Seq_REM = load('all_REM_sweeps_iter_parts.mat', 'all_sweeps');
Seq_REM = Seq_REM.all_sweeps;
Seq_wake = load('all_wake_seqs_parts.mat');
Seq_wake = Seq_wake.Seq_wake;
%% align behav is needed.

[u_tvec, u_idx] = unique(behav.time);
pos(:,1) = interp1(u_tvec,behav.position(u_idx,1),ms.time);
pos(:,2) = interp1(u_tvec,behav.position(u_idx,2),ms.time);
velocity = interp1(u_tvec, behav.speed(u_idx), ms.time);
move_idx = velocity > 2.5;
%% plot an example for track data.v
type = 'parts';
% actual data
this_type = v2_manifold.v2_track(move_idx,:);
this_type = this_type(1:size(Seq_wake.(type).H,2),:); % cut to just the first 3/4 to match the SeqNMF train/test.
data_to_plot = [nan(20,size(Seq_wake.data,2)); Seq_wake.data]; % just for plotting.  fills in with some nans
tvec = ms.time(move_idx)/1000;
tvec = tvec(1:size(Seq_wake.(type).H,2)); % corrected tvec for Seq 3/4


% get seq indices.
if sum(Seq_wake.(type).is_significant) ~=0
    fprintf('Sig seq detected in SeqNMF.  Using window = %is\n', Seq_wake.(type).L);
end
Seq_idx = Seq_wake.(type).H >= prctile(Seq_wake.(type).H, 99.5);

Seq_samples = Seq_wake.(type).L*mode(diff(ms.time)); % convert length of factor into samples.

Seq_idx_out = zeros(size(Seq_idx));
for ii = length(Seq_idx):-1:1
    if Seq_idx(ii)
        Seq_idx_out(ii:ii+Seq_samples) = 1;
    end
end

% get the SEq times
% Seq_per = prctile(Seq_wake.(type).H, 99);

%or try only the top x H values
[H_pks, H_pks_idx] = findpeaks(Seq_wake.(type).H, 'MinPeakDistance', Seq_samples);
[H_top, top_idx] = sort(H_pks, 'descend');
H_top_rel = H_top/max(H_top);
H_pks_idx_sort = H_pks_idx(top_idx(1:5));

Seq_idx_out = zeros(size(Seq_wake.(type).H));
Color_idx = Seq_idx_out;
for ii = 1:length(H_pks_idx_sort)
    Seq_idx_out(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = 1;
    Color_idx(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = ii;
end
Seq_idx_out = logical(Seq_idx_out);
these_idx = find(Seq_idx_out == 1);


% plot everything
figure(1)
subplot(4,8,1)
this_data = squeeze(Seq_wake.(type).W(Seq_wake.(type).indSort,:,:)); % W factor data
imagesc((0:size(Seq_wake.(type).W,3)-1)./ mode(diff(ms.time)),-20:size(Seq_wake.(type).W,1), squeeze(Seq_wake.(type).W(Seq_wake.(type).indSort,:,:)).*([1:size(Seq_wake.(type).W,1)]+(floor(size(Seq_wake.(type).W,1)/2)))');
rectangle('position', [0, -20,(size(Seq_wake.(type).W,3)-1)./ mode(diff(ms.time)),20],'facecolor', 'w', 'edgecolor', [1 1 1])
c_ord = parula(size(H_pks_idx_sort,2));
ylim([-20 size(Seq_wake.(type).W,1)])
text(Seq_wake.(type).L*.5, -10, [upper(type) ' Seq Fac'], 'fontweight', 'bold', 'HorizontalAlignment', 'center')

subplot(4,8,2:8)
% hold on
% plot(Seq_idx, 'r');
% plot(Seq_idx_out, 'b');
hold on
% imagesc(data_to_plot)
imagesc((0:length(Seq_wake.data)-1),1:size(Seq_wake.data,1), Seq_wake.data(Seq_wake.(type).indSort,:).*([1:size(Seq_wake.data,1)]+(floor(size(Seq_wake.data,1)/2)))')
ylabel('cell ID (sorted)')
set(gca,'YDir','reverse')
% axis xy
xlim([0 length(Seq_wake.data)]);
% box off
c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
% make a white block to put seq blocks on top of.
rectangle('position', [0, -20,length(Seq_wake.data)-1,20],'facecolor', 'w', 'edgecolor', [1 1 1])


for ii = 1:length(H_pks_idx_sort)
    rectangle('position', [H_pks_idx_sort(ii), -20,Seq_wake.(type).L*mode(diff(ms.time)),20],'facecolor', c_ord_l(ii,:), 'edgecolor', [1 1 1])
end
% set(gca, 'color', c_ord(1,:))
ylim([-20 size(Seq_wake.data,1)])
% subplot(4,4,[5 9 13])
% title('Sig factor')
% ylabel('cell ID (sorted)')
% imagesc(squeeze(Seq_wake.(type).W(Seq_wake.(type).indSort,:)).*([1:size(Seq_wake.(type).W,1)]+(floor(size(Seq_wake.data,1)/2)))')

subplot(4,8,[9:12, 17:20, 25:28])
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

subplot(4,8, [13:16, 21:24, 29:32])
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
imagesc(1:size(Seq_wake.data,2),1:size(Seq_wake.data,1), Seq_wake.data(Seq_wake.(type).indSort,:).*([1:size(Seq_wake.data,1)]+(floor(size(Seq_wake.data,1)/2)))')
ylabel('cell ID (sorted)')
set(gca,'YDir','reverse')
% axis xy
xlim([min(tvec) max(tvec)]);
axis off
% xlabel('time (s)')

set(gcf, 'position', [440   185   560   660],'MenuBar','none','ToolBar','none','resize','off')
clear F
for ii = 4:2:(33*60*.2)+4
    % for ii =3600:2:5282
    
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
for ii=4:2:(30*60*.2)+4
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(ii) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% manifold movie with Seq using post_REM
% actual data
type = 'parts';
% get seq indices.
sig_thresh = .25;
good_L = zeros(length(Seq_REM.(type)),1);
good_L_val = good_L;
for iSeq = 1:length(Seq_REM.(type))
    if sum(sum(Seq_REM.(type){iSeq}.sig,2))/size(Seq_REM.(type){iSeq}.sig,1) > sig_thresh
        fprintf('Sig seqs (%d/%d) detected in SeqNMF.  Using window = %is\n',sum(sum(Seq_REM.(type){iSeq}.sig,2)),size(Seq_REM.(type){iSeq}.sig,1),Seq_REM.(type){iSeq}.L);
        good_L(iSeq,1) = iSeq;
        good_L_val(iSeq,1) = sum(sum(Seq_REM.(type){iSeq}.sig,2))/size(Seq_REM.(type){iSeq}.sig,1);
    end
end


% get the best L
if sum(good_L,'all') == 0
    error(['No sig factors present at this sig interation threshold (' num2str(sig_thresh) ')']);
end
[good_L_val, sidx] = sort(good_L_val,'descend');

% loop over sig Ls
for iL = 1:sum(good_L ~=0)
    
    this_L = good_L(sidx(iL));
    fprintf('<strong>Using L</strong> = %ds...\n', Seq_REM.(type){this_L}.L)
    
    % get the best iteration (terrible method here).
    sig_idx = 1:length(sum(Seq_REM.(type){this_L}.sig,2)');
    sig_idx = sig_idx(logical(sum(Seq_REM.(type){this_L}.sig,2))');
    [~, best_iter] = max(Seq_REM.(type){this_L}.power(logical(sum(Seq_REM.(type){this_L}.sig,2))'));
    best_iter = sig_idx(best_iter);
    
    
    
    % Prepare raw data using the sig L value
    this_type = v2_manifold.v2_postREM;
    this_type = this_type(1:length(Seq_REM.(type){this_L}.Train),:);
    
    Seq_idx = Seq_REM.(type){this_L}.H{best_iter} >= prctile(Seq_REM.(type){this_L}.H{best_iter}, 99.5);
    
    Seq_samples = Seq_REM.(type){this_L}.L*mode(diff(ms.time)); % convert length of factor into samples.
    
    Seq_idx_out = zeros(size(Seq_idx));
    for ii = length(Seq_idx):-1:1
        if Seq_idx(ii)
            Seq_idx_out(ii:ii+Seq_samples) = 1;
        end
    end
    
    % get the SEq times
    % or try only the top x H values
    [H_pks, H_pks_idx] = findpeaks(Seq_REM.(type){this_L}.H{best_iter}, 'MinPeakDistance', Seq_samples);
    [H_top, top_idx] = sort(H_pks, 'descend');
    H_top_rel = H_top/max(H_top);
    H_pks_idx_sort = H_pks_idx(top_idx(1:5));
    
    Seq_idx_out = zeros(size(Seq_REM.(type){this_L}.H));
    Color_idx = Seq_idx_out;
    for ii = 1:length(H_pks_idx_sort)
        Seq_idx_out(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = 1;
        Color_idx(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = ii;
    end
    Seq_idx_out = logical(Seq_idx_out);
    these_idx = find(Seq_idx_out == 1);
    
    
    % plot raw data with H blocks
    figure(100+iL)
    % info + factor
    subplot(4,8,1)
    this_data = squeeze(Seq_REM.(type){this_L}.W{best_iter}(Seq_REM.(type){this_L}.indSort(:,best_iter),:,:)); % W factor data
    imagesc((0:size(Seq_REM.(type){this_L}.W{best_iter},3)-1)./ mode(diff(ms.time)),-20:size(Seq_REM.(type){this_L}.W{best_iter},1), this_data.*([1:size(this_data,1)]+(floor(size(this_data,1)/2)))');
    rectangle('position', [0, -20,(size(Seq_REM.(type){this_L}.W{best_iter},3)-1)./ mode(diff(ms.time)),20],'facecolor', 'w', 'edgecolor', [1 1 1])
    c_ord = parula(size(H_pks_idx_sort,2));
    ylim([-20 size(Seq_REM.(type){this_L}.W{best_iter},1)])
    text(Seq_REM.(type){this_L}.L*.5, -10, [upper(type) ' Seq Fac'], 'fontweight', 'bold', 'HorizontalAlignment', 'center')
    
    
    
    subplot(4,8,2:8)
    % hold on
    imagesc((0:length(Seq_REM.(type){this_L}.Train)-1)/mode(diff(ms.time)),1:size(Seq_REM.(type){this_L}.Train,1), Seq_REM.(type){this_L}.Train(Seq_REM.(type){this_L}.indSort(:,best_iter),:).*([1:size(Seq_REM.(type){this_L}.Train,1)]+(floor(size(Seq_REM.(type){this_L}.Train,1)/2)))')
    ylabel('cell ID (sorted)');
    
    % set the color profile
    c_ord = parula(size(H_pks_idx_sort,2));
    c_ord_l = linspecer(size(H_pks_idx_sort,2));
    
    
    % plot boxes for the top 5 H values
    for ii = 1:length(H_pks_idx_sort)
        rectangle('position', [H_pks_idx_sort(ii)/mode(diff(ms.time)), -20,Seq_REM.(type){this_L}.L,20],'facecolor', c_ord_l(ii,:), 'edgecolor', [1 1 1])
    end
    set(gca, 'color', c_ord(1,:))
    ylim([-20 size(Seq_REM.(type){this_L}.Train,1)])
    box off
    
    % add vertical lines to separate REM blocks.
    hv = vline(Seq_REM.parts{this_L}.seg_idx/mode(diff(ms.time)));
    for ii = 1:length(hv)
        hv(ii).LineWidth = 2;
        hv(ii).Color = [.9 .9 .9];
    end
    
    subplot(4,8,[9:12, 17:20, 25:28])
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
    
    
    subplot(4,8, [13:16, 21:24, 29:32])
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
    
    mkdir([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session])
    saveas(gcf, [inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'Seq_postREM_manifold_L' num2str(Seq_REM.(type){this_L}.L) '_Iter' num2str(best_iter) '.fig'])
    saveas(gcf,[inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'Seq_postREM_manifold_L' num2str(Seq_REM.(type){this_L}.L) '_Iter' num2str(best_iter) '.png'])
    
end % sig L loop
%% REM movie
% this_lfp = REM_LFP(1:length(this_type)); % limit to training set from Seq data.
% this_lfp = this_lfp./max(this_lfp);

% only use the best L
this_L = good_L(sidx(2));


% get the best iteration (terrible method here).
sig_idx = 1:length(sum(Seq_REM.(type){this_L}.sig,2)');
sig_idx = sig_idx(logical(sum(Seq_REM.(type){this_L}.sig,2))');
[~, best_iter] = max(Seq_REM.(type){this_L}.power(logical(sum(Seq_REM.(type){this_L}.sig,2))'));
best_iter = sig_idx(best_iter);

[H_pks, H_pks_idx] = findpeaks(Seq_REM.(type){this_L}.H{best_iter}, 'MinPeakDistance', Seq_samples);
[H_top, top_idx] = sort(H_pks, 'descend');
H_top_rel = H_top/max(H_top);
H_pks_idx_sort = H_pks_idx(top_idx(1:5));

Seq_idx_out = zeros(size(Seq_REM.(type){this_L}.H));
Color_idx = Seq_idx_out;
for ii = 1:length(H_pks_idx_sort)
    Seq_idx_out(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = 1;
    Color_idx(H_pks_idx_sort(ii):H_pks_idx_sort(ii)+Seq_samples) = ii;
end
Seq_idx_out = logical(Seq_idx_out);
these_idx = find(Seq_idx_out == 1);



% set up the movie plot
c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
tvec = 1:length(this_type);
tvec = tvec*(1/mode(diff(ms.time)));

this_data = Seq_REM.(type){this_L}.Train;

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


subplot(8,1, 6:8)
imagesc((0:length(this_data)-1)/mode(diff(ms.time)),1:size(this_data,1), this_data(Seq_REM.(type){this_L}.indSort(:,best_iter),:).*([1:size(this_data,1)]+(floor(size(this_data,1)/2)))')
ylabel('cell ID (sorted)')
set(gca,'YDir','reverse')
% axis xy
xlim([0 (length(this_data)-1)/mode(diff(ms.time))]);
% xlabel('time (s)')

c_ord = parula(size(H_pks_idx_sort,2));
c_ord_l = linspecer(size(H_pks_idx_sort,2));
for ii = 1:length(H_pks_idx_sort)
    rectangle('position', [H_pks_idx_sort(ii)/mode(diff(ms.time)), -20,Seq_REM.(type){this_L}.L,20],'facecolor', c_ord_l(ii,:), 'edgecolor', [1 1 1])
end
set(gca, 'color', c_ord(1,:))
ylim([-20 size(Seq_REM.(type){this_L}.Train,1)])
box off
% add vertical lines to separate REM blocks.
hv = vline(Seq_REM.parts{this_L}.seg_idx/mode(diff(ms.time)));
for ii = 1:length(hv)
    hv(ii).LineWidth = 2;
    hv(ii).Color = [.9 .9 .9];
end


% subplot(8,1,6)
% imagesc(tvec, 1, nan(size(tvec)))
% axis off
% y_lim = ylim;
% xlabel('time (s)')


clear F
% for ii = 4:2:length(Seq_idx_out)
for ii = 8000:1:length(Seq_idx_out)
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
    
    
    %     subplot(8,1,6)
    %     if ii< 34
    %         xlim([tvec(1) tvec(66)])
    %     else
    %         xlim([tvec(ii)-1 tvec(ii)+1])
    %     end
    %     ylim(y_lim);
    % add time.
    %     t1 = text(10, .5, [num2str(tvec(ii),3) 's']);
    %
    %     xlim([min(pos(:,1)) max(pos(:,1))]);
    %     ylim([-1 1]); set(gca, 'ytick', []);
    %     xlabel('position (cm)')
    %     set(gca, 'XDir','reverse')
    subplot(8,1,6:8)
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
writerObj = VideoWriter([inter_dir filesep 'Seq_manifold_EC' filesep subject filesep session filesep 'REM_full_L' num2str(Seq_REM.(type){this_L}.L) '_Iter' num2str(best_iter) '.avi']);
writerObj.FrameRate = round(16);
writerObj.Quality = 90;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
% for ii=4:2:(30*60*2)+4
% for ii = 4:2:length(F)
for ii = 8000:1:length(Seq_idx_out)

    % convert the image to a frame
    frame = F(ii) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


