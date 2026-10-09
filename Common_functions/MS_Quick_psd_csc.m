function psd = MS_Quick_psd_csc(csc, Chan_to_use)
%% MS_Quick_psd_csc: generates PSDs for each csc channel 
%   
% EC 20236-10-08 initial verstion
%   Based on Quick_psd by EC back in the van der Meer lab. 


if nargin < 2
    Chan_to_use = 1:lenth(csc.length); 
end


csc = MS_select_tsd(csc, Chan_to_use); 
%% set up defaults. 

line_width = 2;

c_ord = linspecer(length(csc.label)); % predefine some colours.

%% get the psd
cfg_psd = [];
cfg_psd.hann_win = 2^12; % always make this in base 2 for speed

for iSite = 1:length(csc.label)
    [psd.(csc.label{iSite}).pxx, psd.(csc.label{iSite}).f] = pwelch(csc.data(iSite,:), hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , csc.cfg.hdr{1}.SamplingFrequency);
end


%% plot the PSD for each channel on one plot.  

color.blue = double([158,202,225])/255;
color.green = double([168,221,181])/255;
color.red = [0.9    0.3    0.3]; 

subplot(3,2,[2,3,5,6])
% sort channel labels
for iSite = 1:length(csc.label)
    this_num = regexp(csc.label{iSite},'\d*','Match');
    num_labels(iSite) = str2num(this_num{1});
end
[~, sorted_labels_idx] = sort(num_labels); 

l = 0; 
for iSite = sorted_labels_idx
    l = l+1;
    hold on
    plot3(psd.(csc.label{iSite}).f, 10*log10(psd.(csc.label{iSite}).pxx),ones(size(psd.(csc.label{iSite}).f))*iSite, 'color', c_ord(l,:),'linewidth', line_width);
end
xlim([0 250])
y_val = ylim;
xlabel('Frequency (Hz)')
% colour bars for specific frequencies of interest. 
legend(strrep(csc.label(sorted_labels_idx), '_', ' '),  'location', 'southwest', 'orientation', 'vertical');
rectangle('position', [1, y_val(1), 4, y_val(2) - y_val(1)],  'facecolor', [color.red 0.2], 'edgecolor', [color.red 0.2])
rectangle('position', [6, y_val(1), 4, y_val(2) - y_val(1)],  'facecolor', [color.blue 0.2], 'edgecolor', [color.blue 0.2])
rectangle('position', [40, y_val(1), 30, y_val(2) - y_val(1)],  'facecolor', [color.green 0.2], 'edgecolor', [color.green 0.2])
%flip obj order for clearity
chi=get(gca, 'Children');
set(gca, 'Children',flipud(chi))

if isunix
d_name = strsplit(cd, '/'); % what to name the figure.  
else
    d_name = strsplit(cd, '\'); % what to name the figure.  
end
title(strrep(d_name{end}, '_', ' ')) 


% zoom in on 0-14 hz in a small plot
axes('Position',[.7 .6 .2 .3])
box on
% add colours for freq bands
l = 0; 
for iSite = sorted_labels_idx
    hold on
    l = l+1; 
    plot(psd.(csc.label{iSite}).f, 10*log10(psd.(csc.label{iSite}).pxx), 'color', c_ord(l,:),'linewidth', line_width);
end
xlim([0 12])
y_val = ylim;

rectangle('position', [1, y_val(1), 4, y_val(2) - y_val(1)],  'facecolor', [color.red 0.2], 'edgecolor', [color.red 0.2])
rectangle('position', [6, y_val(1), 4, y_val(2) - y_val(1)],  'facecolor', [color.blue 0.2], 'edgecolor', [color.blue 0.2])

set(gca,'yticklabels', [])
% flip obj order
chi=get(gca, 'Children');
set(gca, 'Children',flipud(chi))
pos = get(gcf, 'position');
set(gcf, 'position', [pos(1) pos(2) pos(3)*1.5 pos(4)*2])
% %% plot the coherence between pairs of Chan_to_use. 
% subplot(2,3,[4,5, 6])
% for iPairs = 1:length(labels)
%     hold on
%     plot(coh.(labels{iPairs}).f, coh.(labels{iPairs}).p, 'color', c_ord(length(Chan_to_use)+iPairs,:), 'linewidth', line_width)
% end
% xlim([0 120])
% ylim([0 1])
% legend(pair_labels)%, 'location', 'EastOutside', 'orientation', 'vertical')
% maximizesfn 2019sfn 2019
% saveas(gcf,[save_dir filesep 'PSD_check.png'])
