%% shifts test


d_mat = [];

for ii = 1:length(ms.shifts)
    
    this_d = cell2mat({ms.shifts{ii}(:).diff});
    
    d_mat = [d_mat this_d]; 
end


% not helpful

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


plot(flash_enc, h.Values(flash_t < ms.time(ms.timestamps(1))), 'xb')
plot(flash_rec, h.Values(flash_t > ms.time(sum(ms.timestamps(1:2))+2)), 'xr')

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

%%
% move_idx 

mov_idx = behav_rec_a.speed > 2.5; 


figure(910)
clf
hold on
plot(behav_rec.time , behav_rec.position(:,1))
plot(ms_rec.time , ms_rec.denoise(:,1))

MS_event_rate_map(flash_vec(mov_idx) , ms_rec.time(mov_idx), behav_rec_a.position(mov_idx,:), 1, 1);





