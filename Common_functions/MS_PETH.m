function [peth] = MS_PETH(cfg_in, S, t)
%% KA_PETH: creates a simple histogram of spikes across some time intervals.
%
%
%  Follows methods from Fraser et al. 2023 (https://www.biorxiv.org/content/10.1101/2023.06.28.546936v1.full)
%
% 'Peri-stimulus time histograms (PSTHs) were constructed around
% event-related responses using 0.01 ms bins. The spiking activity of each
% neuron across these bins of the PSTH was smoothed using a half-normal
% filter (σ = 6.6) that used activity in previous, but not upcoming, bins.
% To visualize the normalized activity of neurons, the mean activity within
% each of the smoothed bins of the PSTH was transformed to a z-score as
% follows: (Fi – Fmean)/ FSD, where Fi is the firing rate of the ith bin of
% the PSTH, and Fmean and FSD are the mean and SD of the firing rate of the
% 10 s baseline period. Color-coded maps and average traces of individual
% neurons’ activity were constructed based on these z-scores.'
%
%
%  Inputs:
%     -cfg_in: [struct]   contains user configuration parameters. will
%     override default fields.
%
%     - S: [struct]    spike data in the TS format. Will only process one
%     cell in index 1.
%
%     - IV: [struct]  contains intervals in the IV format.
%
%
%  Outputs:
%     -peth: [struct] contains the histogram of spikes along with the cfgs.
%



%% initialize

cfg_def = [];
cfg_def.bin = 0.01; % assuming seconds not '0.01ms' as per Fraiser.
cfg_def.b_line = -10; % 10s baseline.
cfg_def.sd = 6.6; % as per Fraiser.
cfg_def.gauss_window = 1;
cfg_def.gauss_sd = 0.02; 

cfg_def.plot = 0; % toggle for plots.
cfg_def.evt_color_mat = repmat([0 0 0], length(t),1);
cfg_def.markersize = 5;


cfg = ProcessConfig(cfg_def, cfg_in);

%% convert times (t) to IV

% IV = iv(t + cfg.win(1), t+cfg.win(2));

%% Align spikes to trials.



% set up gau kernal
gauss_window = cfg.gauss_window./cfg.bin; % 1 second window
gauss_SD = cfg.gauss_sd./cfg.bin; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg.bin; % normalize by binsize



S_out = []; outputT = []; outputS = [];
for iT = length(t):-1:1
    S0 = restrict(S, t(iT)+cfg.win(1), t(iT)+cfg.win(2));


    tbin_edges = t(iT)+cfg.win(1):cfg.bin:t(iT)+cfg.win(2);
    tbin_centers = tbin_edges(1:end-1)+cfg.bin/2;


    outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
    outputS = [outputS; S0.t{1}-t(iT)];

        S_gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
        if size(S_gau_sdf,1) >1
            S_gau_sdf = S_gau_sdf';
        end
        outputGau(:,iT) = S_gau_sdf;



    % get the histogram around the events.

    [spk_count, B] = histcounts(S0.t{1},tbin_edges);
    spk_count =  spk_count(1:end-1);

    S_out(iT,:) = spk_count; % save the counts.
end


S_rate = sum(S_out,1)./(tbin_edges(end) - tbin_edges(1));




mean_S_gau = nanmean(outputGau,2); % get the mean gaussian smoothed firing rate

mean_S_gau_z = nanmean(outputGau,2); % get the mean gaussian smoothed firing rate
mean_S_gau_z = (mean_S_gau_z - nanmean(outputGau_shuf,2))./nanstd(outputGau_shuf, [], 2);

idx = nearest_idx3(0, outputIT); % get the event time index

pre_stim_means = nanmean(outputGau(1:idx-1,:),1);  % get the mean of the gau smoothed firing rate before the event.
post_stim_means = nanmean(outputGau(idx:end,:),1); % % get the mean of the gau smoothed FR after the event.

pre_stim_std = nanstd(outputGau(1:idx-1,:),[],1);  % get the mean of the gau smoothed firing rate before the event.
post_stim_std = nanstd(outputGau(idx:end,:),[],1); % % get the mean of the gau smoothed FR after the event.

%% apply filtering


peth = [];
%% plotting if needed.

subplot(2,1,1);
% 	imagesc(window,[1 nT], outputID);
% 	colormap(1-0.25*gray);
hold on;
u_val = unique(outputT);
for iV = 1:length(u_val)
    this_idx = outputT == u_val(iV);
    if isempty(this_idx)
        continue
    end
    plot(outputS(this_idx), outputT(this_idx)+0.5,'.', 'color', cfg.evt_color_mat(u_val(iV),:), 'MarkerSize', cfg.markersize)
end

%     plot(outputS, outputT+0.5, '.k', 'MarkerSize', 5);
xlabel('peri-event (sec)');
ylabel('Event #');
if length(t) < 2
    ylim([0 1])
else
    ylim([1 length(t)])
end



    subplot(2,1,2);
    yyaxis left
    hold on
plot(outputIT, mean_S_gau,'color', 'k', 'linewidth', cfg.linewidth)
        xlim(cfg.window);
            ylabel('firing rate (Hz)');
         if ~isempty(cfg.shuff)
            plot(outputIT,nanmean(outputGau_shuf,2), '--', 'color', [0.3 .3 .3])
        end
        
        if ~(max(mean_S_gau)) ==0
            ylim([min(mean_S_gau) max(mean_S_gau)])
            if (size(t,2) > 1) || (size(t,1) == 1)
                rectangle('position', [0 min(mean_S_gau) 0.001  abs(max(mean_S_gau))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            else
                rectangle('position', [0 min(mean_S_gau) abs(mode(t(:,2)-t(:,1)))  abs(max(mean_S_gau))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            end
        end
    end
