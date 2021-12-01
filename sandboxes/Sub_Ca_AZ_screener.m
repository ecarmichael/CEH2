% temp sub checker


%% plot the traces 

cfg = [];
cfg.Ca_chan = 1:15%[1 2 4 5 7 13 15 16 17 18 19 21 22 23 24]; % get a subset of cells.
cfg.plot_type = '2d';
cfg.offset = 4;
MS_plot_ca(cfg, data.ms)

xlim([0 600])

%% check some sig cells
p_vals = [];
for ii = length(data.SI):-1:1
p_vals(ii) = data.SI(ii).spatial.place.MI_pval;
end

[p_sort_val, p_sort] = sort(p_vals, 'ascend');

% 
% for ii = 1:20
%     figure(ii)
%     imagesc(data.SI(ii).spatial.place.posterior)
%     title(num2str(ii))
% end
    
%% nice plots
for iC = 136 %p_sort(1:25)
%     move_idx = (data.behav_aligned.speed < 25) & (data.behav_aligned.speed > 2); 
move_idx =1:length(data.behav_aligned.position(:,1));
spk_idx = find(data.ms.Binary(move_idx,iC));
figure(iC*100)
subplot(1,3,1)
hold on
plot(data.behav_aligned.position(move_idx,1),data.behav_aligned.position(move_idx,2), 'color', [.7 .7 .7])
plot(data.behav_aligned.position(spk_idx,1),data.behav_aligned.position(spk_idx,2), 'r.')
xlim([min(data.behav_aligned.position(move_idx,1)), max(data.behav_aligned.position(move_idx,1))]);
ylim([min(data.behav_aligned.position(move_idx,2)), max(data.behav_aligned.position(move_idx,2))])
title(iC)
subplot(1,3,2)
  imagesc(data.SI(iC).spatial.place.occupany')
set(gca, 'ydir', 'normal', 'xdir', 'normal');
title(num2str(iC))

subplot(1,3,3)
  imagesc(imgaussfilt(data.SI(iC).spatial.place.posterior', .5))
set(gca, 'ydir', 'normal', 'xdir', 'normal');

Square_subplots
end

