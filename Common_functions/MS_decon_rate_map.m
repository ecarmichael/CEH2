function [rate_map, occ_mat] = MS_decon_rate_map(data,tvec,pos,X_bins,Y_bins, plot_flag, smooth_SD)
%% MS_decon_rate_map:  estimates the firing rate map form the deconvolved Ca signal. Requires the OASIS toolbox
%
%
%
%    Inputs: 
%    - data: [1 x nSamples]  deconvolved signal (from OASIS).   
%
%    - tvec: [1 x nSamples]  time vector. 
%
%    - pos [1 or 2 x nSamples]  position data as x;y
%
%    - bin_size: bin size in cm. 
%
%    - plot_flag [binary]  if 1 plot the outputs.  default is 0. 
%
%    - smooth_SD  number of pixels for the SD.  typically 2
%
%
%
%
%    Outputs: 
%    - rate_map  [Xbins x Ybins] rate map
%
%    - occ_map  [Xbins x Ybins] occupancy map
%
%
%
%
% EC 2021-05-14   initial version 
%   2021-05-15  - Switched to imgaussfilt instead of kernal + conv2
%
%
%% initialize

if nargin <3
   error('Requires data, tvec and position data')
   
elseif nargin < 4
X_bins = 0:bin_size:ceil(max(pos(:,1)));
% X_bin_centers = X_bins +  bin_size/2;
% X_bin_centers = X_bin_centers(1:end-1);

% same for Y bins
Y_bins = 0:bin_size:ceil(max(pos(:,2)));
% Y_bin_centers = Y_bins +  bin_size/2;
plot_flag = 0; 
    
elseif nargin < 6
    plot_flag = 0; 
end

if nargin <7
%       smooth_bin = [0, 0]; 
      smooth_SD = 0; 
end

% %% split position into bins

% Y_bin_centers = Y_bin_centers(1:end-1);


% add in some guassian smoothing. 
% kernel = gausskernel(smooth_bin,smooth_SD); % Gaussian kernel of 4x4 pixels, SD of 2 pixels (note this should sum to 1)

%% plot for sanity

if plot_flag
    spk_idx = data > 0; 
    if ishandle(101)
        close
    end
    figure(101)
    subplot(2,2,1)
    hold on
    plot(pos(:,1), pos(:,2), 'color', [.8 .8 .8])
    plot(pos(spk_idx,1), pos(spk_idx,2), '.r', 'markersize',8)
end


% get occupancy. 

% get n frames of occupancy per bin. 
if size(pos,2) == 2
    occ_mat= zeros(length(X_bins)-1,length(Y_bins)-1);
    spk_mat= zeros(length(X_bins)-1,length(Y_bins)-1);
else
    occ_mat = zeros(length(X_bins)-1);
end

% get occupancy and event rate per bin 
for iY = 1:length(Y_bins)-1
    for iX = 1:length(X_bins)-1
        this_vec = zeros(size(data)); 
        p_idx = find(pos(:,1) >= X_bins(iX) & pos(:,1) < X_bins(iX+1) & pos(:,2) >= Y_bins(iY) & pos(:,2) < Y_bins(iY+1));
        
        if ~isempty(p_idx)
        this_data = data(p_idx); % get all the data points while in this bin
        spk_mat(iX, iY) = sum(this_data > 0); 
        occ_mat(iX, iY) = length(p_idx)/length(tvec); 
        end
    end
end

% smooth if need (will use zeros if not specified on input)
% occ_mat = conv2(occ_mat, kernel, 'same'); 
% spk_mat = conv2(spk_mat, kernel, 'same');
 
occ_mat = imgaussfilt(occ_mat, smooth_SD); 
spk_mat = imgaussfilt(spk_mat, smooth_SD); 
        

no_occ_idx = find(occ_mat == 0); 
occ_mat(no_occ_idx) = NaN; % Nan out bins that were not visited.
spk_mat(no_occ_idx) = NaN; % same for 

occ_mat = occ_mat./(tvec(end)/1000); % convert to time in s from tvec. 

rate_map = spk_mat./occ_mat; % compute tuning curve. 


% plot the occupancy matrix
if plot_flag
    subplot(2,2,2)
    imagesc(X_bins,Y_bins,  occ_mat');
    set(gca, 'YDir', 'normal', 'xDir', 'normal');
    title('occupancy')
    h = colorbar;
    ylabel(h, 'secs');
end   


% plot the 'spike' matrix 
if plot_flag
    subplot(2,2,3)
    imagesc(X_bins, Y_bins, spk_mat');
        set(gca, 'YDir', 'normal', 'xDir', 'normal');
    title('spikes')
    h = colorbar;
    ylabel(h, 'nspikes');
end

% plot the tuning curve
if plot_flag
    subplot(2,2,4)
    imagesc(X_bins, Y_bins, rate_map');
        set(gca, 'YDir', 'normal', 'xDir', 'normal');
    title('rate map')
    h = colorbar;
    ylabel(h, 'spk/s');
end




