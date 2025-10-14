%% JC_REM_Fig1_LFP_EMG


data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1';
nlx_dir  = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Jisoo\Raw_LFP\2019-07-08_09-03-55_PV1069_LTD1'; 
cd(data_dir)


%% load the LFP data

load([data_dir 'ms_resize.mat'])

cd(nlx_dir)

cfg.fc = {'CSC1.ncs' 'CSC7.ncs'};
cfg.desired_sampling_frequency = 1000; 
csc = MS_LoadCSC(cfg);

%% cut out some REM epochs
rem_idx = find(contains(ms_seg_resize.hypno_label, 'REM')); 

csc_rem = restrict(csc, ms_seg_resize.NLX_csc{rem_idx(end)}.tvec(1)-60, ms_seg_resize.NLX_csc{rem_idx(end)}.tvec(end)+60);

Fs = round(1/mode(diff(csc_rem.tvec)));

c_ord = MS_linspecer(3); 
tvec = csc_rem.tvec-csc_rem.tvec(1); 
%% basic plot
figure(101)
clf

ax(1) = subplot(5,1,1);
plot(tvec, csc_rem.data(2,:), 'k')
xline([60 tvec(end)-60], 'color', c_ord(2,:))

ax(2) = subplot(5,1,2);
plot(tvec, csc_rem.data(1,:),'k')
xline([60 tvec(end)-60], 'color', c_ord(2,:))


ax(3) = subplot(5,1,3:5);
fb = cwtfilterbank(SignalLength=numel(csc_rem.data(2,:)),SamplingFrequency=Fs,...
    FrequencyLimits=[1 120], VoicesPerOctave=32);
cwt(csc_rem.data(2,:), FilterBank=fb);
[cfs,frq] = cwt(csc_rem.data(2,:), FilterBank=fb);
% [cfs,frq] = cwt(csc_rem.data(2,:), Fs, 'bump', 'VoicesPerOctave',32);
% helperCWTTimeFreqPlot(cfs,tvec,frq, 'surf','CWT of Bat Echolocation','seconds','Hz')

cfs_pow = abs(cfs).^2;
[minf,maxf] = cwtfreqbounds(numel(csc_rem.data(2,:)),Fs);
freq = 2.^(round(log2(minf)):round(log2(maxf)));

% imagesc(tvec, 10*log10(frq), cfs_pow)

% xlabels = get(gca, 'xtick');
% set(gca, 'xticklabels', xlabels - win_s)
h = gcf;
x_lim = xlim;
hold on
% xline(win_s, '--w', 'start', 'linewidth', 2);
%         xline(x_lim(2) - win_s, '--k', 'end', 'linewidth', 2);
yline(10, '--w', '10hz', 'linewidth', 2, 'LabelHorizontalAlignment','left');

AX = gca;

AX.YTickLabelMode = 'auto';
AX.YTick = freq;
ylim([1 140])


%%
% dt is the sampling period in unit of time and T is the length of your recording
T = 100;
dt = 0.001;
t = 0:dt:T;
NumVoices = 64;
a0 = 2^(1/NumVoices);
wavCenterFreq = centfrq('morl');

minfreq = 0.5;
maxfreq = 120;

minscale = wavCenterFreq/(maxfreq*dt);
maxscale = wavCenterFreq/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales = a0.^(minscale:maxscale).*dt;

% testing the wavelet spectrogram on a sinasoid with frequency of 150 Hz
% cwt_test = cwtft({sin(2*pi*150.*t),dt},'scales',scales,'wavelet','morl');
[cwt_test, f] = cwt(csc_rem.data(2,:),'amor', Fs, VoicesPerOctave=32, FrequencyLimits=[minfreq maxfreq]);
cwt_test_power = abs(cwt_test).^2;
% cwt_test_power = abs(cwt_test.cfs).^2;
% freq = cwt_test.frequencies;
% freq = f;
figure;
% args = {t,freq,10*log10(cwt_test_power)};
% surf(args{:},'edgecolor','none');
% view(0,90);
imagesc(tvec, f, cwt_test_power);
set(gca, 'YDir', 'normal')
axis tight;
colormap('jet');

title('wavelet spectrogram');
xlabel('time (s)');
ylabel('frequency (Hz)');
h = colorbar;
set(h, 'YTick', get(gca, 'CLim'));
ylabel(h, 'power (db)', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'Bold', ...
    'Position', [1.6304 -23.5000 0]);