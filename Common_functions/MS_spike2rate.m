function [rate_tsd] = MS_spike2rate(S, tvec, bin_s, sig)




if nargin < 2
    error('Needs a spike TS and a tvec')
elseif nargin < 3
    bin_s = 0.001; 
    sig = 0.05;
elseif nargin < 4
    sig = 0.05;
end


rate_tsd = tsd; 

% construct internal tvec
dt = mode(diff(tvec));
t_edges = tvec(1):dt:tvec(end);
t_centers = t_edges(1:end-1)+dt/2;

%% loop over Spikes

for ii = length(S.t):-1:1

    spk_binned = histc(S.t{ii},t_edges);
    spk_binned(end-1) = spk_binned(end-1)+spk_binned(end); % combine last bin and edge
    spk_binned = spk_binned(1:end-1);

    % construct gaussian kernel
    gauss_window = 1./dt;
    gauss_SD = sig./dt; % 0.02 seconds (20ms) SD
    gk = gausskernel(gauss_window,gauss_SD); gk = gk./dt; % normalize by binsize

    % convolve
    rate_tsd.data(ii,:) = conv2(spk_binned,gk,'same');
    
    rate_tsd.label = S.label; 
end

rate_tsd.tvec = t_centers'; 
rate_tsd.units = 'Hz'; 
rate_tsd.cfg.history.mfun{end+1,1} = 'MS_spike2rate'; 
rate_tsd.cfg.history.cfg{end+1,1}.dt = dt; 
rate_tsd.cfg.history.cfg{end,1}.sig = sig; 
rate_tsd.cfg.history.cfg{end,1}.bin_s = bin_s; 
