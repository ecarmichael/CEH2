function [psd, coh] = MS_Quick_psd(Chan_to_use, type)
%% MS_Quick_psd: generates PSDs for each .ncs in the current directory specified in 'Chan_to_use' 
%   (if empty do all of them). Uses pwelch method.use 'type' input as either 'long' to use all of 
%   the data or 'fast' to use the first 1/4. 
%
%
%   Inputs:
%       - Chan_to_use: [1xnCells] Cells containing the file names. Eg:
%       Chan_to_use = {'CSC1.ncs', 'CSC2.ncs'}. 
%       - type: [string] can be either 'long' (all data) or 'fast' (first
%       1/4) of the data). Default: 'long'
%
%
%   Outputs:
%       - psd: [struct] power sepctral densities and frequencies from pwelch for each
%       specified channel. along with cfg field with parameters.  
%
% EC 2020-06-01 initial verstion
%   Based on Quick_psd by EC back in the van der Meer lab. 


%% set up defaults. 

if nargin == 0
    Chan_to_use = {};
    type = 'long';
    disp('''Chan_to_use'' not specified.  Using all .ncs files')
    disp('''type'' no specified, using ''long'' for all data. Use type = ''fast'' for only the first 1/4 of the data to save time')

elseif nargin == 1
    type = 'long';
    disp(' ''type'' not specified, using ''long'' for all data. Use type = ''fast'' for only the first 1/4 of the data to save time')
end

line_width = 2;
%% Loop through channels specified in the ExpKeys and get the PSD.  
% Chan_to_use = Chan_to_use;
cfg = [];
cfg.fc = Chan_to_use; 
% cfg.decimateByFactor = 2; % leave empty for reactive decimation using cfg.desired_sampling_frequency
cfg.desired_sampling_frequency = 2000; % helps with speed. 

csc = MS_LoadCSC(cfg);

%% get the psd
cfg_psd = [];
cfg_psd.hann_win = 2^12; % always make this in base 2 for speed
switch type
    case 'fast'
        cut = round(length(csc.data(1,:))/4); % this can be used for speed.
        %     only takes in half the data.
        for iSite = 1:length(csc.label)
            [psd.(csc.label{iSite}(1:end-4)).pxx, psd.(csc.label{iSite}(1:end-4)).f] = pwelch(csc.data(iSite,1:cut), hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , csc.cfg.hdr{1}.SamplingFrequency);
        end
        
    case 'long'
        for iSite = 1:length(csc.label)
            [psd.(csc.label{iSite}(1:end-4)).pxx, psd.(csc.label{iSite}(1:end-4)).f] = pwelch(csc.data(iSite,:), hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*2 , csc.cfg.hdr{1}.SamplingFrequency);
        end
        
end
%
% %% Get the coherence between pairs of Chan_to_use. [move this to it's own
% function!!!] 
% if length(Chan_to_use) >1
%     site_comb = nchoosek(Chan_to_use, 2); % get all combinations of Chan_to_use
%     labels = [];
%     for iComb = 1:length(site_comb); % loop through combinations of Chan_to_use.  
%         label = [site_comb{iComb, 1}(1:end-4) '_' site_comb{iComb, 2}(1:end-4)];
%         S1 = strfind(Chan_to_use, site_comb{iComb,1}); % find the corresponding
%         S2 = strfind(Chan_to_use, site_comb{iComb,2}); % find the other corrsponding site idx
%         S1 = find(not(cellfun('isempty', S1))); % gets the actual idx
%         S2 = find(not(cellfun('isempty', S2)));
%         % get teh coherence using the same paramters as in pwelch.  
%         [coh.(label).p, coh.(label).f] =mscohere(csc.data(S1,:), csc.data(S2,:),hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*4 , csc.cfg.hdr{1}.SamplingFrequency);
%         labels{iComb} = label;
%         pair_labels{iComb} = [real_label{S1}(1:3) '-' real_label{S2}(1:3)] % get the name of the pair.  
%     end
% end

%     %%
%
%
% for iSite = 1:length(Chan_to_use)-1;
%     label = [Chan_to_use{iSite}(1:end-4) '_' Chan_to_use{iSite+1}(1:end-4)];
%     [coh.(label).p, coh.(label).f] =mscohere(csc.data(iSite,1:cut), csc.data(iSite+1,1:cut),hanning(cfg_psd.hann_win), cfg_psd.hann_win/2, cfg_psd.hann_win*4 , csc.cfg.hdr{1}.SamplingFrequency);
%     labels{iSite} = label;
%     pair_labels{iSite} = [real_label{iSite}(1:3) '-' real_label{iSite+1}(1:3)]
% end
%% make a plot of the targets [BETA] relative to bregma. 

% electrode plot (beta)
if isfield(ExpKeys, 'Target')
    figure
c_ord = linspecer(length(csc.label)); % predefine some colours. 

    subplot(2,3,1)
    targets = fieldnames(ExpKeys.Target);
    hold on
    for iT = 1:length(targets)
        flag = 0;
        for iLab = 1:length(ExpKeys.Chan_to_use_labels)
            if strcmp(ExpKeys.Chan_to_use_labels{iLab}(1:3), targets{iT})
                flag = 1;
            end
        end
        if flag ==1
            rectangle('Position',[ExpKeys.Target.(targets{iT})(2)-.1, ExpKeys.Target.(targets{iT})(1)-.25, .5, .5],'Curvature',[1 1], 'facecolor', [c_ord(iT,:), .3])
        else
            rectangle('Position',[ExpKeys.Target.(targets{iT})(2)-.1, ExpKeys.Target.(targets{iT})(1)-.25, .5, .5],'Curvature',[1 1])
        end
        text(ExpKeys.Target.(targets{iT})(2), ExpKeys.Target.(targets{iT})(1),-ExpKeys.Depths.(targets{iT})(1), targets{iT})
    end
    text(0, 0, 'Bregma')
    xlim([-2 4])
    ylim([-2 2])
    zlim([-10 0])
end

%% plot the position data.note fails if there is more than one ..nvt file (can happen in RR2)
% subplot(2,3,4)
% if exist('VT1.zip') % unzips the nvt if it is zipped.  
%     unzip('VT1.zip')
% end
% cfg_pos = [];
% pos = LoadPos(cfg_pos);
% plot(pos.data(1,:), pos.data(2,:), '.')
% axis off
%% plot the PSD for each channel on one plot.  
subplot(3,2,[2,3,5,6])
for iSite = 1:length(csc.label)
    hold on
    plot(psd.(csc.label{iSite}(1:end-4)).f, 10*log10(psd.(csc.label{iSite}(1:end-4)).pxx), 'color', c_ord(iSite,:),'linewidth', line_width);
end
xlim([0 200])
y_val = ylim;
xlabel('Frequency (Hz)')
% colour bars for specific frequencies of interest. 
% theta (broad 4-12)
% rectangle('position', [4, y_val(1), 8, y_val(2)-y_val(1)], 'facecolor', [0 0.2 .8 0.2],'edgecolor', [0 0.2 .8 0.2]); %low gamma rectangle
% low gamma
% rectangle('position', [45, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [0 0.8 .2 0.2],'edgecolor', [0 0.8 .2 0.2]) %high gamma rectangle
% high gamma
% rectangle('position', [70, y_val(1), 20, y_val(2)-y_val(1)], 'facecolor', [0 0.8 .2 0.2],'edgecolor', [0 0.8 .2 0.2]) %high gamma rectangle

legend(strrep(csc.label, '_', ' '),  'location', 'southwest', 'orientation', 'vertical');


if isunix
d_name = strsplit(cd, '/'); % what to name the figure.  
else
    d_name = strsplit(cd, '\'); % what to name the figure.  
end
title(strrep(d_name{end}, '_', ' ')) 


% zoom in on 0-14 hz in a small plot
axes('Position',[.7 .6 .2 .3])
box on
for iSite = 2:length(csc.label)
    hold on
    plot(psd.(csc.label{iSite}(1:end-4)).f, 10*log10(psd.(csc.label{iSite}(1:end-4)).pxx), 'color', c_ord(iSite,:),'linewidth', line_width);
end
xlim([0 14])
y_val = ylim;
set(gca,'yticklabels', [])


% %% plot the coherence between pairs of Chan_to_use. 
% subplot(2,3,[4,5, 6])
% for iPairs = 1:length(labels)
%     hold on
%     plot(coh.(labels{iPairs}).f, coh.(labels{iPairs}).p, 'color', c_ord(length(Chan_to_use)+iPairs,:), 'linewidth', line_width)
% end
% xlim([0 120])
% ylim([0 1])
% legend(pair_labels)%, 'location', 'EastOutside', 'orientation', 'vertical')
% maximize
saveas(gcf, 'PSD_check.png')
