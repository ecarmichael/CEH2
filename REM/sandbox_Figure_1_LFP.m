%% JC_REM_Fig1_LFP_EMG


data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1';
nlx_dir  = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Jisoo\Raw_LFP\2019-07-08_09-03-55_PV1069_LTD1'; 
save_dir = 'C:\Users\ecarm\Desktop\REM_figs';
cd(data_dir)


%% load the LFP data

load([data_dir filesep 'ms_resize.mat'])

cd(nlx_dir)

cfg.fc = {'CSC1.ncs' 'CSC7.ncs'};
cfg.desired_sampling_frequency =4000; 
csc = MS_LoadCSC(cfg);

%% cut out some REM epochs
rem_idx = find(contains(ms_seg_resize.hypno_label, 'REM')); 
win = 60; 
this_idx = 7;
csc_rem = restrict(csc, ms_seg_resize.NLX_csc{rem_idx(this_idx)}.tvec(1)-win, ms_seg_resize.NLX_csc{rem_idx(this_idx)}.tvec(end)+win);

Fs = round(1/mode(diff(csc_rem.tvec)));

c_ord = MS_linspecer(3); 
tvec = csc_rem.tvec-csc_rem.tvec(1); 

%% spectrogram

win_size = 2^11; % works out to 512 samples. Out of some superstition I always use base 2 (b/c of bytes or something) when computing spectra. 
n_overlap = win_size/4; % just needs to be smaller than the window. 1/4 gives nice temporal resolution. 
fs = csc.cfg.hdr{1}.SamplingFrequency;
freq_range = 1:0.1:30; % range of frequencies. 

[~,F,T, P] = spectrogram(csc_rem.data(2,:),win_size,n_overlap,freq_range,fs); 

% T = T+csc_rem.tvec(1); % just so that our time vectors are the same
% for ii = size(P,2):-1:1
    P_s = smoothdata(P,2, 'gaussian',5); 
% end
%% basic plot
z_win = 2; % zoom window
figure(101)
clf

ax(1) = subplot(6,1,1);
plot(tvec, csc_rem.data(2,:), 'k')
xline([win tvec(end)-win], 'color', c_ord(2,:))
ylim([-0.001 0.001])
axis off

ax(2) = subplot(6,1,2);
plot(tvec, csc_rem.data(1,:),'k')
xline([win tvec(end)-win], 'color', c_ord(2,:))
ylim([-0.001 0.001])
axis off

ax(3) = subplot(6,1,3:4);
cla
args = {T,F,10*log10(P_s)};
surf(args{:},'edgecolor','none');
view([0 90])
ylim([F(1) F(end)])
set(gca, 'YDir', 'normal')
hold on
    rectangle('position', [win+10, F(1), z_win, F(end)-F(1)], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth',2)
    rectangle('position', [win-40, F(1), z_win, F(end)-F(1)], 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',2)
    rectangle('position', [floor(tvec(end)-win+5), F(1), z_win, F(end)-F(1)], 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',2)

linkaxes(ax, 'x')
xlim([win-60 tvec(end)-win+60])
SetFigure([], gcf, 1)


%% samples
figure(102)
clf

% NREM zoom
ar(1) = subplot(6,3,1);
plot(tvec, csc_rem.data(2,:), 'k')
ylim([-0.001 0.001]); axis off

ar(2) = subplot(6,3,4);
plot(tvec, csc_rem.data(1,:),'k')
ylim([-0.001 0.001]); axis off

ar(3) = subplot(6,3,[7 10]);
cla
args = {T,F,10*log10(P_s)};
surf(args{:},'edgecolor','none');
view([0 90])
set(gca, 'YDir', 'normal')
ylim([F(1) F(end)]); box off

linkaxes(ar, 'x')
xlim([win-40 win-40+z_win])

%REM zoom
ay(1) = subplot(6,3,2);
plot(tvec, csc_rem.data(2,:), 'k')
ylim([-0.001 0.001]); axis off

ay(2) = subplot(6,3,5);
plot(tvec, csc_rem.data(1,:),'k')
ylim([-0.001 0.001]); axis off

ay(3) = subplot(6,3,[8 11]);
cla
args = {T,F,10*log10(P_s)};
surf(args{:},'edgecolor','none');
view([0 90])
set(gca, 'YDir', 'normal')
ylim([F(1) F(end)])

linkaxes(ay, 'x')
xlim([win+10 win+10+z_win])

%wake zoom
ay(1) = subplot(6,3,3);
plot(tvec, csc_rem.data(2,:), 'k')
ylim([-0.001 0.001]); axis off

ay(2) = subplot(6,3,6);
plot(tvec, csc_rem.data(1,:),'k')
ylim([-0.001 0.001]); axis off

ay(3) = subplot(6,3,[9 12]);
cla
args = {T,F,10*log10(P_s)};
surf(args{:},'edgecolor','none');
view([0 90])
set(gca, 'YDir', 'normal')
ylim([F(1) F(end)]); box off

linkaxes(ay, 'x')
xlim([floor(tvec(end)-win+5) floor(tvec(end)-win+5)+z_win])

SetFigure([], gcf, 1)

%% save the figures
figure(101)
print("-bestfit",[save_dir filesep strrep(ms_seg_resize.dirName(4:end), '\', '_') '_fig1_LFP'], '-dpdf')

figure(102)
print("-bestfit",[save_dir filesep strrep(ms_seg_resize.dirName(4:end), '\', '_') '_fig1_LFP_inset'], '-dpdf')

%%



fb = cwtfilterbank(SignalLength=numel(csc_rem.data(2,:)),SamplingFrequency=Fs,...
    FrequencyLimits=[1 120], VoicesPerOctave=32);
cwt(csc_rem.data(2,:), FilterBank=fb);
[cfs,frq] = cwt(csc_rem.data(2,:), FilterBank=fb);
% [cfs,frq] = cwt(csc_rem.data(2,:), Fs, 'bump', 'VoicesPerOctave',32);
% helperCWTTimeFreqPlot(cfs,tvec,frq, 'surf','CWT of Bat Echolocation','seconds','Hz')

cfs_pow = abs(cfs).^2;
[minf,maxf] = cwtfreqbounds(numel(csc_rem.data(2,:)),Fs);
freq = 2.^(round(log2(minf)):round(log2(maxf)));


figure(103)
bx(3) = subplot(5,1,3:5);
% imagesc(tvec, frq, 10*log10(cfs_pow))
args = {tvec,frq,10*log10(cfs_pow)};
surf(args{:},'edgecolor','none');
view(0,90);
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
% ylim([1 140])


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
args = {tvec,f,10*log10(cwt_test_power)};
surf(args{:},'edgecolor','none');
view(0,90);
% imagesc(tvec, f, cwt_test_power);
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

%% 

minfreq = 0.5;
maxfreq = 30;

 cwt(csc_rem.data(2,:),'amor', Fs, VoicesPerOctave=32, FrequencyLimits=[minfreq maxfreq]);
% [cfs, frq] = cwt(csc_rem.data(2,:), fs); 
%  xlabels = get(gca, 'xtick'); 
%  set(gca, 'xticks', xlabels/60)
        % 
         AX = gca;
        % [minf,maxf] = cwtfreqbounds(numel(csc_rem.data(2,:)),fs);
        % 
        % freq = 2.^(round(log2(minf)):round(log2(maxf)));
        % AX.YTickLabelMode = 'auto';
        % AX.YTick = freq;
% ylim([1 140])
% Color gradient
numColors = 256;
colormap('parula');
xlim([])