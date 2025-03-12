%% shifts test


d_mat = [];

for ii = 1:length(ms.shifts)
    
    this_d = cell2mat({ms.shifts{ii}(:).diff});
    
    d_mat = [d_mat this_d]; 
end


% not helpful

%% try getting the median of only the negative values by removing all positives before getting the median. 


temp_detrend = NaN(size(ms.detrendRaw));
p_mat = []; 
for ii = size(ms.detrendRaw,2):-1:1
    
    temp_detrend(:,ii) = ms.detrendRaw(:,ii);
    
    
    [~, p_idx] = findpeaks(-temp_detrend(:,ii), 'MinPeakHeight', .2, 'MinPeakDistance', 10);
    if length(p_idx) > 100
    
        p_mat = [p_mat, p_idx']; 
        keep_idx(ii) = true; 
    temp_detrend(temp_detrend(:,ii) > 0,ii) = 0; 
    
    else
        keep_idx(ii) = false; 
    end
    
    
    
end

%%
figure(161)

clf
ax(1) = subplot(2,1,1);
plot(ms.time, median(temp_detrend(:,keep_idx),2))

ax(2) = subplot(2,1,2);
hold on
h = histogram(ms.time(p_mat), round(length(ms.time)/5));
yline(prctile(h.Values, 99.5))
linkaxes(ax, 'x')

% get the indicies of the putative flashes from the histogram. 


% flash_idx = 

%% grab specific frames. 