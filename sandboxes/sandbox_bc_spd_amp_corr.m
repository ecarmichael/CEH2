%% sandbox_BC_spd_amp_corr


% Load some data
cfg_fc = [];
cfg_fc.fc = {'CSC2.ncs'}; 
csc = MS_LoadCSC(cfg_fc);

pos = MS_DLC2TSD(cd, [], [7.2 7.2])

% filter and convert to amplitude
 cfg_filt_t = [];
    cfg_filt_t.type = 'cheby1';                                            %'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt_t.f  = [4 12];                                                % freq range to match Mizuseki et al. 2011
    cfg_filt_t.order = 3;                                                  %type filter order
    cfg_filt_t.display_filter = 0;                                         % use this to see the fvtool
    theta_csc = FilterLFP(cfg_filt_t, csc);                                % filter the raw LFP using
    theta_csc.data = abs(hilbert(theta_csc.data));    

    
%% decimate to match sampling frequency
cfg_deci = [];
cfg_deci.decimateFactor = 20; 

theta_amp_d = decimate_tsd(cfg_deci,theta_csc); 


% interpolate the spd to match the downsampled LFP amp

spd_interp = interp1(pos.tvec, pos.data(5,:), theta_amp_d.tvec, 'next', 'extrap'); 

%% check the interp/downsample


figure(91919)
clf
hold on
plot(csc.tvec, csc.data(1,:))
plot(theta_csc.tvec, theta_csc.data(1,:), 'r')
plot(theta_amp_d.tvec, theta_amp_d.data(1,:), '--m')
% 
yyaxis right
plot(pos.tvec, pos.data(5,:), 'k')
plot(theta_amp_d.tvec, spd_interp, '--g')

%% compute the xcorr between spd_interp and theta_amp_d

nan_idx = isnan(spd_interp) | isnan(theta_amp_d.data'); % find any NaN values in either data

spd_t_amp_corr = corr(spd_interp(~nan_idx), theta_amp_d.data(~nan_idx)'); 



%% histogram  of spd x amplitude


