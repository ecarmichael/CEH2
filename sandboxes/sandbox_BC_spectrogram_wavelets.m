%% sandbox_BC_spec_wav


%% load some data

cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_LINEAR_TRACK\raw\BC1807_07_12_2023_D4_INHIBITION');

cfg_csc = [];
cfg_csc.fc = {'CSC4.ncs'};
csc = MS_LoadCSC(cfg_csc);

cfg_csc = [];
cfg_csc.fc = {'CSC1.ncs'};
emg = MS_LoadCSC(cfg_csc);


evts = LoadEvents([]);

laser_on = evts.t{4};
laser_off = evts.t{3}; % assuming one is on and the other is off. 
%% plot data to see if the evts look right
figure(101)
clf
subplot(3,1,1)
cla
hold on


for ii = 1:length(laser_on)
    rectangle('position', [laser_on(ii), min(csc.data), laser_off(ii) - laser_on(ii), max(csc.data) - min(csc.data)], 'FaceColor', 'g', 'EdgeColor', 'g')

end

plot(csc.tvec, csc.data)

xlim([csc.tvec(1) csc.tvec(end)])

%% make one big spectrogram
win_size = 2^9; % works out to 512 samples. Out of some superstition I always use base 2 (b/c of bytes or something) when computing spectra. 
n_overlap = win_size/4; % just needs to be smaller than the window. 1/4 gives nice temporal resolution. 
fs = csc.cfg.hdr{1}.SamplingFrequency;
freq_range = 1:0.25:120; % range of frequencies. 

[~,F,T, P] = spectrogram(csc.data,win_size,n_overlap,freq_range,fs); 

T = T+csc.tvec(1); % just so that our time vectors are the same

subplot(3,1,2)
cla
imagesc(T, F, 10*log10(P))
set(gca, 'YDir', 'normal')
hold on
for ii = 1:length(laser_on)
    rectangle('position', [laser_on(ii), 0, laser_off(ii) - laser_on(ii), 120], 'FaceColor', 'none', 'EdgeColor', 'g')

end


%% same thing but wavelet;

% I made a wrapper function for the built in cwt wavelet function. This
% will let you put the EMG and raw LFP on top. but you can also open it up
% and pull parts from it. 


[cfs, frq] = cwt(csc.data, fs); 
%  xlabels = get(gca, 'xtick'); 
%  set(gca, 'xticks', xlabels/60)

         AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(csc.data),fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
ylim([1 140])
%% wavelet for the whole thing. 

figure
[cfs, frq]=cwt(csc_r.data, fs);
     AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(csc_r.data),fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
ylim([1 140])
%% event sepcific

for ii = 1:length(laser_on)
    
    this_idx = nearest_idx(laser_on(ii), csc.tvec); 
    figure(ii)
MS_plot_cw_emg(csc.tvec, csc.data, emg.data, this_idx, 5);
xline(laser_on(ii) - laser_off(ii))

end

%% possible CWT diff
csc_on = restrict(csc, laser_on(1), laser_off(1)); 
csc_off = restrict(csc, laser_off(1), laser_off(1) - laser_on(1)); 

 figure
[cfs, frq]=cwt(csc_on.data, fs);
     AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(csc_on.data),fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
ylim([1 140])


%% get an event triggered spectrogram showing the mean around the laster offset

% this requires the fieldtrip toolbox: https://github.com/fieldtrip/fieldtrip
% once you have it you can add it to the path and run ft_defaults.  Its a
% good idea to remove it from your path using 'rmpath(<path to ft>)'
% because it overwrite some functions like 'butter'

% this is a wrapper for the FT function. The inputs are just the CSC, some
% event times, a label, a frequency of interest range, baseline to
% normalize tom and then a timewindows. I set the window to be big for the
% wavelet cone and then restrict it a bit closer with an xlim. 
Triggered_Spec_FT(csc, laser_on, 'laser on', 1:.1:100, [-2 0], [-5 5])
xlim([-2 2])
