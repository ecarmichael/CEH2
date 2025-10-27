%% Figure 1 LFP and Ca raw trace example for a sleep expisode. 


if ismac
usr_dir = '/Users/ecar/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter'; 
else
usr_dir= 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter'; 

end
data_dir = strrep([usr_dir filesep 'Sample_data_code\Figure 1'], '\', filesep);
nlx_dir  = strrep([usr_dir filesep  'Sample_data_code\Figure 1\nlx_data'], '\', filesep); 
save_dir = strrep([usr_dir filesep  'Desktop\REM_figs'], '\', filesep);
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
this_idx = 1;
csc_rem = restrict(csc, ms_seg_resize.NLX_csc{rem_idx(this_idx)}.tvec(1)-win, ms_seg_resize.NLX_csc{rem_idx(this_idx)}.tvec(end)+win);

Fs = round(1/mode(diff(csc_rem.tvec)));

c_ord = MS_linspecer(3); 
tvec = csc_rem.tvec-csc_rem.tvec(1); 

%% spectrogram

win_size = 2^11; % works out to 512 samples. Out of some superstition I always use base 2 (b/c of bytes or something) when computing spectra. 
n_overlap = win_size/2; % just needs to be smaller than the window. 1/4 gives nice temporal resolution. 
fs = csc.cfg.hdr{1}.SamplingFrequency;
freq_range = 1:0.1:30; % range of frequencies. 

[~,F,T, P] = spectrogram(csc_rem.data(2,:),win_size,n_overlap,freq_range,fs); 

% smooth for visualization

P_s = smoothdata(P,2, 'gaussian',5); 
%% basic plot
z_win = 1; % zoom window
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
    rectangle('position', [win+25, F(1)-1, z_win, F(end)-F(1)+1], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth',2)
    rectangle('position', [win-30, F(1)-1, z_win, F(end)-F(1)+1], 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',2)
    rectangle('position', [floor(tvec(end)-win+5), F(1)-1, z_win, F(end)-F(1)+1], 'FaceColor', 'none', 'EdgeColor', 'k','LineWidth',2)

cb = colorbar; 
c_pos = cb.Position; 
cb.Position = c_pos+[.05 0 0 0]; 
cb.Ticks = [cb.Ticks(1) cb.Ticks(end)]; 

linkaxes(ax, 'x')
xlim([win-60 tvec(end)-win+60])
SetFigure([], gcf, 1)
colormap(magma)


%% samples
figure(102)
clf

% NREM zoom
ar(1) = subplot(6,3,1);
plot(tvec, csc_rem.data(2,:), 'k')
ylim([-0.0005 0.0005]); axis off

ar(2) = subplot(6,3,4);
plot(tvec, csc_rem.data(1,:),'k')
ylim([-0.0005 0.0005]); axis off

ar(3) = subplot(6,3,[7 10]);
cla
args = {T,F,10*log10(P_s)};
surf(args{:},'edgecolor','none');
view([0 90])
set(gca, 'YDir', 'normal')
ylim([F(1) F(end)]); box off

linkaxes(ar, 'x')
xlim([win-30 win-30+z_win])

%REM zoom
ay(1) = subplot(6,3,2);
plot(tvec, csc_rem.data(2,:), 'k')
ylim([-0.0005 0.0005]); axis off

ay(2) = subplot(6,3,5);
plot(tvec, csc_rem.data(1,:),'k')
ylim([-0.0005 0.0005]); axis off

ay(3) = subplot(6,3,[8 11]);
cla
args = {T,F,10*log10(P_s)};
surf(args{:},'edgecolor','none');
view([0 90])
set(gca, 'YDir', 'normal')
ylim([F(1) F(end)])

linkaxes(ay, 'x')
xlim([win+25 win+25+z_win])

%wake zoom
ay(1) = subplot(6,3,3);
plot(tvec, csc_rem.data(2,:), 'k')
ylim([-0.0005 0.0005]); axis off

ay(2) = subplot(6,3,6);
plot(tvec, csc_rem.data(1,:),'k')
ylim([-0.0005 0.0005]); axis off

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
colormap(magma)

%% add a plot for the calium

figure(103)
ax(3) = subplot(3,2,[3 5]);
cla

data = ms_seg_resize.RawTraces{this_idx}(); 
time = ms_seg_resize.time{this_idx}; 
cfg.offset =.5; 
cfg.width = 1; 
c_ord = viridis(96); 
cfg.rescale =  'a'; 
tick_val = []; 
hold on
for iC = 1:64
    if strcmp(cfg.rescale, 'zscore')
        plot(time*0.001, zscore(data(:,iC))+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
        % plot3(time*0.001, repmat(iC,size(data,1),1), zscore(data(:,iC))+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
        % tick_val(end+1) = median(zscore(data(:,iC))+iC*cfg.offset);
    elseif strcmp(cfg.rescale, 'max')
        plot(time*0.001, data(:,iC)./max(data(:,iC))+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
        % tick_val(iC) = median(data(:,iC)./max(data(:,iC))+iC*cfg.offset);
    else
        plot(time*0.001, data(:,iC)+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
        % tick_val(iC) = median(data(:,iC)+iC*cfg.offset);
    end
    % tick_label{iC} = iC;
end
    
axis off
ylim([min(tick_val)-1 max(tick_val)+cfg.offset*2])
xlim([25 35])
SetFigure([], gcf, 1)

%% save the figures
figure(101)
print("-bestfit",[save_dir filesep strrep(ms_seg_resize.dirName(4:end), '\', '_') '_fig1_LFP'], '-dpdf', "-vector")

figure(102)
print("-bestfit",[save_dir filesep strrep(ms_seg_resize.dirName(4:end), '\', '_') '_fig1_LFP_inset'], '-dpdf', "-vector")

figure(103)
print("-bestfit",[save_dir filesep strrep(ms_seg_resize.dirName(4:end), '\', '_') '_fig1_REM_ca_inset'], '-dpdf', "-vector")

    