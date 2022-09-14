%% sandbox_HPP_Test

data_dir = '/media/williamslab/Baldr/SWR_detector_test_1252_maze3/2022-09-13_SWR_tests';


cd(data_dir)


%% load data

cfg_lfp = [];
cfg_lfp.fc = {'CSC4.ncs'};

csc = MS_LoadCSC(cfg_lfp);


evts = LoadEvents([]);

rec_dur = [];
for ii = length(evts.t{1}):-1:1
   rec_dur(ii) = evts.t{2}(ii) - evts.t{1}(ii);
    
    
end


csc_r = restrict(csc, evts.t{1}(end-1), evts.t{2}(end-1)); 
evts_r = restrict(evts, evts.t{1}(end-1), evts.t{2}(end-1)); 


pulse = sort([evts_r.t{3} evts_r.t{4} evts_r.t{5}]); 
pulseIV = iv(pulse, pulse+0.025); 
%% detect SWRs in test recording

cfg_swr = [];
cfg_swr.check = 0; % plot checks.
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [150 250]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter =0; % use this to see the fvtool


% smoothing
cfg_swr.kernel.samples = csc_r.cfg.hdr{1}.SamplingFrequency/100;
cfg_swr.kernel.sd = csc_r.cfg.hdr{1}.SamplingFrequency/100;

% detection
cfg_swr.threshold =1.5;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
cfg_swr.merge_thr = 0.01; %merge events that are within 20ms of each other.
%         cfg_swr.nan_idx = SWD_idx; % where are any nans, say from excluding artifacts, other events...

% restrictions
cfg_swr.max_len = [];
cfg_swr.max_len.operation = '<';
cfg_swr.max_len.threshold = .1;

%                 cfg_swr.min_len = [];
%                 cfg_swr.min_len.operation = '<';
%                 cfg_swr.min_len.threshold = .2;
cfg_swr.nCycles = 20; % number of cycles
cfg_swr.nCycles_operation = '<='; % number of cycles

% variaence
cfg_swr.var = [];
cfg_swr.var.operation = '<';
cfg_swr.var.threshold = 1;

[SWR_evts, SWR_filt, SWR_amp] = MS_get_LFP_events_sandbox(cfg_swr, csc_r);

cfg_plot.display = 'tsd';
PlotTSDfromIV(cfg_plot, SWR_evts, csc_r)

SWR_centers = IVcenters(SWR_evts); 

%% limit to event that overlap 

hit_IV = IntersectIV([], SWR_evts, pulseIV);

hit_centers = IVcenters(hit_IV); 

fprintf('Hit Rate %.2f%% (%d/%d)\n', (length(hit_centers)/length(SWR_evts.tstart))*100, length(hit_centers), length(SWR_evts.tstart)) 
%% generate a histogram of pulse times relative to SWR onset. 

for ii = length(pulseIV.tstart):-1:1
    
    this_range = -(SWR_evts.tend - pulseIV.tstart(ii)); % find the event that occured closest to the end fo the SWR. 
    keep = this_range <=0; 
    min_t = min(-this_range(keep)); 
    this_iv = find(this_range == -min_t); 
    
    off_t(ii) = pulseIV.tstart(ii) - SWR_evts.tstart(this_iv); 
    
end

off_t(off_t <0) = []; 

figure(101)
histogram(off_t*1000, 25, 'BinLimits' , [0 100])
xlim([0 100])
xline(mean(off_t)*1000, 'linewidth', 5, 'label', ['mean = ' num2str(mean(off_t)*1000,'%.2f\n') 'ms']); 
xlabel('time from SWR onset (ms)')
ylabel('Count (hits only)')
SetFigure([], gcf)
%% plot the average LFP at time of pulse
win_s = .5;
Fs = csc_r.cfg.hdr{1}.SamplingFrequency;

swr_filt = []; swr_amp = []; 

for ii = length(hit_centers):-1:1
    this_idx = nearest_idx3(hit_centers(ii), csc_r.tvec); 
    if this_idx+(win_s*Fs) > length(csc_r.tvec) || this_idx-(win_s*Fs) <1
        continue
    else
    swr_filt(ii,:) = SWR_filt.data(1,this_idx-(win_s*Fs):this_idx+(win_s*Fs));   
    swr_amp(ii,:) = SWR_amp.data(1,this_idx-(win_s*Fs):this_idx+(win_s*Fs));   

    end
    
end

figure(1)
clf
hold on
yyaxis left
plot(-.5:(1/Fs):.5, nanmean(swr_filt))
ylabel('mean ripple band (150-250Hz)')
yyaxis right
plot(-.5:(1/Fs):.5, nanmean(swr_amp));
ylabel('mean ripple amplitude')

xline(0, '--k', 'Pulse', 'linewidth', 3)
xlabel('Time from SWR candidate event detection (s)')
set(gca, 'xtick', -.5:.1:.5)
SetFigure([], gcf)
maximize