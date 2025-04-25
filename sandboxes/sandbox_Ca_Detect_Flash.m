%% shifts test

% 
% d_mat = [];
% 
% for ii = 1:length(ms.shifts)
%     
%     this_d = cell2mat({ms.shifts{ii}(:).diff});
%     
%     d_mat = [d_mat this_d]; 
% end
% 
% 
% % not helpful

%% try getting the median of only the negative values by removing all positives before getting the median. 
ms = MS_Ca_good_cells(ms); 

temp_detrend = NaN(size(ms.detrendRaw));
p_mat = []; keep_idx= [];
for ii = size(ms.detrendRaw,2):-1:1
    
    temp_detrend(:,ii) = diff([ms.detrendRaw(:,ii);  ms.detrendRaw(end,ii)]);
    
    
    [~, p_idx] = findpeaks(-temp_detrend(:,ii), 'MinPeakHeight', .1, 'MinPeakDistance', 10);
    if length(p_idx) > 100
    
        p_mat = [p_mat, p_idx']; 
        keep_idx(ii) = true; 
    temp_detrend(temp_detrend(:,ii) > 0,ii) = 0; 
    
    else
        keep_idx(ii) = false; 
    end
    
    
    
end

keep_idx = logical(keep_idx);
%%
figure(161)

med_detrend = median(temp_detrend(:,keep_idx),2); 

clf
ax(1) = subplot(2,1,1);
hold on
plot(ms.time, med_detrend)
% plot(ms.time(sort(p_mat)), med_detrend(sort(p_mat)), 'x')


ax(2) = subplot(2,1,2);
hold on
h = histogram(ms.time(p_mat), round(length(ms.time)/5));
yline(prctile(h.Values, 99))

% get the indicies of the putative flashes from the histogram. 

bin_centers = h.BinEdges(1:end-1)+h.BinWidth; 

flash_idx =h.Values > prctile(h.Values, 99);
flash_t = bin_centers(flash_idx);

flash_enc = flash_t(flash_t < ms.time(ms.timestamps(1)));
flash_rec = flash_t(flash_t > ms.time(sum(ms.timestamps(1:2))+2));


plot(flash_enc-(h.BinWidth/2), h.Values(flash_t < ms.time(ms.timestamps(1))), 'xb')
plot(flash_rec-(h.BinWidth/2), h.Values(flash_t > ms.time(sum(ms.timestamps(1:2))+2)), 'xr')

xline(ms.time(sum(ms.timestamps(1:2))+3))
linkaxes(ax, 'x')



%% grab specific frames. 
% restrict to recall only

ms_rec = MS_restrict(ms, ms.time(sum(ms.timestamps(1:2))+3), flash_t(end));

flash_idx = nearest_idx(flash_rec, ms_rec.time);

flash_vec = zeros(size(ms_rec.time));

flash_vec(flash_idx) = 1;
flash_vec = logical(flash_vec);




flash_rec = flash_rec -ms_rec.time(1);

ms_rec.time = ms_rec.time - ms_rec.time(1); 

behav_rec_a = MS_align_data(behav_rec, ms_rec);
%% hack to get the trials periods
cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\JKA_HPC_07')
Radial_log_JKA_07_2025_02_20_Day5
cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\JKA_HPC_07\2025_02_20')

trl_s = nearest_idx(Rad.D2025_02_20.JKA_07.recall.tstart, behav_rec_a.time); 
trl_e = nearest_idx(Rad.D2025_02_20.JKA_07.recall.tend, behav_rec_a.time); 

trl_idx = zeros(size(behav_rec_a.time)); 
for ii = 1:length(trl_s)
trl_idx(trl_s(ii):trl_e(ii)) = 1; 

end
trl_idx = logical(trl_idx); 

figure(101)
clf
subplot(2,2, 3:4)
hold on
plot(behav_rec_a.time, behav_rec_a.position, 'color', [.7 .7 .7])
plot(behav_rec_a.time(trl_idx), behav_rec_a.position(trl_idx,:))

subplot(2,2,1)
hold on
% plot(behav_rec_a.position(~trl_idx,1), behav_rec_a.position(~trl_idx,2), '.', 'color', [.7 .7 .7]))

plot(behav_rec_a.position(trl_idx,1), behav_rec_a.position(trl_idx,2), '.r')

subplot(2,2,2)
hold on
plot(behav_rec_a.position(~trl_idx,1), behav_rec_a.position(~trl_idx,2), '.', 'color', [.7 .7 .7])
% plot(behav_rec_a.position(trl_idx,1), behav_rec_a.position(trl_idx,2), '.r')



%%
% move_idx 

mov_idx = (behav_rec_a.speed > 2.5) & (behav_rec_a.speed <20);

k_idx = mov_idx & trl_idx; 

figure(910)
clf
hold on
plot(behav_rec.time , behav_rec.position(:,1))
plot(ms_rec.time , ms_rec.denoise(:,1))

MS_event_rate_map(flash_vec(k_idx) , ms_rec.time(k_idx), behav_rec_a.position(k_idx,:), 2, 3);

title(['Rewarded arms: ' num2str(Rad.D2025_02_20.correct)])



%% make a movement movie

this_behav = behav_rec_a; 

figure(888)
clf
plot(this_behav.position(:,1), this_behav.position(:,2), '.', 'color', [.8 .8 .8])

for ii = 1:length(this_behav.time)
    
    
    hold on
    if trl_idx(ii) == 1
             h = plot(this_behav.position(ii,1), this_behav.position(ii,2), 'or', 'markersize', 15);
    else
            h = plot(this_behav.position(ii,1), this_behav.position(ii,2), 'ob', 'markersize', 15);
    end
    
    t = text(80, 80, num2str(this_behav.time(ii)));
    
    drawnow
    pause(mode(diff(this_behav.time))/10)
    delete(h)
    delete(t)
end



