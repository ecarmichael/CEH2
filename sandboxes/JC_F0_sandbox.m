%% dF/F per neuron

load('ms.mat')

%% bin the data using a little hack for non-overlapping bins
Fs = mode(diff(ms.time));
iChan = 25
bins = 1:(Fs*10):length(ms.RawTraces(:,iChan)); 
var_binned = NaN(size(ms.RawTraces(:,iChan))); 
var_vec = NaN(size(ms.RawTraces(:,iChan))); 
for ii = 1:length(bins)
    if ii == length(bins)
        var_binned(ii) = var(ms.RawTraces(bins(ii):end,iChan));
        var_vec(bins(ii):end) = var_binned(ii); 
    else
        var_binned(ii) = var(ms.RawTraces(bins(ii):bins(ii+1),iChan)); 
        var_vec(bins(ii):bins(ii+1)) = var_binned(ii);
    end
end

%% Get the variance and interpolate

p_cut = prctile(var_binned, 20);

% visualize
histogram(var_binned); 
xline(p_cut,'r', '20prct');

F0 = ms.RawTraces(:,iChan); 
F0(var_vec>p_cut) = NaN; 


figure(909)
subplot(4,1,1)
h = plot(ms.time, ms.RawTraces(:,iChan), ms.time, F0);
title('Raw trace')
legend({'Raw Trace', 'Pseudo-baseline'}); 
c = h(2).Color; 
F0 = fillmissing(F0, 'linear'); 

subplot(4,1,2)
plot( ms.time, F0, 'color', c)
title('F0 trace')
legend({'interp F0'}); 

subplot(4,1,3)
plot(ms.time, var_vec); 
legend({'original variance'})
title('variance')

subplot(4,1,4)
plot(ms.time, var_vec); 
ylim([0 p_cut*10])
ylabel('var')
xlabel('time (ms)')
title('zoom in')
yline(p_cut)
legend({'original variance', 'var cut'})





%% %%%%%%%%%% %% dF/F for the entire image

%% grab some raw data and concatenate it. 
vid1= VideoReader('msCam1.avi'); 
data1= squeeze(read(VideoReader('msCam1.avi'))); 
data2= squeeze(read(VideoReader('msCam2.avi'))); 
data3= squeeze(read(VideoReader('msCam3.avi'))); 
data4= squeeze(read(VideoReader('msCam4.avi'))); 
data5= squeeze(read(VideoReader('msCam5.avi'))); 
data6= squeeze(read(VideoReader('msCam6.avi'))); 

dat_cat = double(cat(3, data1, data2, data3, data4, data5, data6)); % convert to double to make var() work. 

clear data*
% infer time from framerate
tvec = 0:1/vid1.FrameRate:((length(dat_cat))/vid1.FrameRate); 
tvec = tvec(1:end-1); % correct for last frame. 


%% bin the data and grab the var

t_bins = 0:10:tvec(end); 
idx_bins = 1:vid1.FrameRate*10:length(dat_cat); 

var_binned = NaN(size(t_bins));
var_vec = NaN(size(tvec)); 
for ii = 1:length(idx_bins) % probably more effient way to do 3d binning, but this works.
    if ii == length(idx_bins)
        var_binned(ii) = var(dat_cat(:,:,idx_bins(ii):end), 0, 'all'); % binned value. 
        var_vec(idx_bins(ii):end) = var_binned(ii);                    % full vector with binned values inserted. Used later. 
    else
        var_binned(ii) = var(dat_cat(:,:,idx_bins(ii):idx_bins(ii+1)), 0, 'all');
        var_vec(idx_bins(ii):idx_bins(ii+1)) = var_binned(ii);
    end
end

%% get the bins in lower 20 percentile
p_cut = prctile(var_binned, 20);

% visualize
histogram(var_binned, length(idx_bins)); 
xline(p_cut,'r', '20prct');

%% interp over bins in upper 80 percentile
var_cut_vec = var_vec; 
var_cut_vec(var_vec>p_cut) = NaN; 

var_cut_vec = fillmissing(var_cut_vec, 'linear'); 

figure(909)
plot(tvec, var_vec, tvec, var_cut_vec); 
legend({'original variance', 'interp var'})


% dat_cat_binned = reshape(dat_cat, vid1.FrameRate, []);


