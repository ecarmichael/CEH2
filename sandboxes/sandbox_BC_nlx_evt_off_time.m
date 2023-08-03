%% sandbox  BC laser pulse start-end



%% load some data
evts = LoadEvents([]); 

cfg =[];
cfg.fc = {'CSC6.ncs'};

csc = MS_LoadCSC(cfg);


%% get some events of interest

on_idx = find(contains(evts.label, '(0x0002)'));  % find the correct label. 

on_ts = evts.t{on_idx}; % grab those times


off_idx = find(contains(evts.label, '(0x0000)'));  % find the corresponding off label
off_ts = evts.t{off_idx}; 

% convert to 'iv' (inverval) format for simplicity

laser_iv = iv(on_ts, off_ts); 


%% alternative. Use a little function based on the cell above. 
pattern = 'TTL Input on AcqSystem1_0 board 0 port 3 value (0x0002).'; 

iv_out = MS_get_evts_off(evts, pattern); 

%% plot some examples
figure(101)
clf

hold on

max_csc = max(csc.data); 
min_csc = min(csc.data); 

offset = (max_csc - min_csc); % this makes the high of the rectangle. 

% c_ord = linspecer(length(on_ts));
Artch_green = [0.530 0.820 0.645];
for ii = 1:length(laser_iv.tstart)
    
    rectangle('position', [laser_iv.tstart(ii), min_csc, laser_iv.tend(ii) - laser_iv.tstart(ii), offset], 'facecolor', [Artch_green .2], 'edgecolor', [Artch_green .2])
    
end

% plot the data on top
plot(csc.tvec, csc.data);


xlim([csc.tvec(1) csc.tvec(end)])


%% gettting the power and phase

% filter the LFP in the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 12]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool

theta_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using

theta_amp = abs(hilbert(theta_csc.data)); % get the amplitude

theta_phi  = angle(hilbert(theta_csc.data(1,:))); 

theta_csc.data = theta_csc.data(1,:); 



%% same but for slow gamma

% filter the LFP in the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [30 58]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 1; % use this to see the fvtool

SG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using

SG_amp = abs(hilbert(SG_csc.data)); % get the amplitude




%% same but for slow gamma

% filter the LFP in the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [70 90]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 4; %type filter order
cfg_filt_t.display_filter = 1; % use this to see the fvtool

FG_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using

FG_amp = abs(hilbert(FG_csc.data)); % get the amplitude


%% Try some Phase-amp coupling


mod_th_g = MS_ModIdx_win(theta_csc, SG_csc, 30*theta_csc.cfg.hdr{1}.SamplingFrequency);
